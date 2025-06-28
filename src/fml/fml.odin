/*
MIT License

Copyright (c) 2025 Rowan Apps

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the “Software”), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/*
FML - Fluid Markup Language

FML is a JSON like configuration format with an added expression parser for host evalauted expressions.

See the README for syntax and usage.

Lots referened from: https://github.com/odin-lang/Odin/blob/master/core/encoding/json/unmarshal.odin
*/

package fml

import sa "core:container/small_array"
import "core:fmt"
import "core:log"
import "core:mem"
import vmem "core:mem/virtual"
import "core:reflect"
import "core:slice"
import "core:strings"

// The strings in any of the errors are not copied, they reference the file data directly.
// so make sure to handle the error, or clone the string, before you free file data.
Unmarshal_Error :: union {
	Unexpected_Char,
	Unexpected_Token,
	Type_Error, //schema and file have differing types for a key
	Missing_Key, //expected a key, didnt find in file
	Duplicate_Key, //key defined > 1 times in the file
	Unknown_Key, //key has no matching entry in schema
	Unsupported_Type, //schema is not allowed to unmarshal to this type
	mem.Allocator_Error,
}

Unexpected_Char :: struct {
	got:      byte,
	expected: Maybe(byte),
	location: Location,
}
Unexpected_Token :: struct {
	got:      Token,
	expected: Maybe(Token_Kind),
}
Type_Error :: struct {
	got:      Token,
	expected: typeid,
}
Missing_Key :: struct {
	key:   string,
	table: string,
}
Duplicate_Key :: struct {
	key:   string,
	table: string,
}
Unknown_Key :: struct {
	key:   string,
	table: string,
}
Unsupported_Type :: struct {
	tid: typeid,
}

// unmarshal `data` into the schema provided by `v`, `path` is for error reporting.
// strings and slices will cause allocations.
// in the case of an error, some of the fields may be allocated, where others have not.
// So check for nil before freeing, or use an arena.
unmarshal :: proc(data: []u8, path: string, v: any, allocator := context.allocator) -> Unmarshal_Error {
	context.allocator = allocator

	v := v
	v = reflect.any_base(v)
	ti := type_info_of(v.id)
	if !reflect.is_pointer(ti) || ti.id == rawptr do return Unsupported_Type{ti.id}

	root_table_obj := any{(^rawptr)(v.data)^, ti.variant.(reflect.Type_Info_Pointer).elem.id}
	p := parser_init(data, path) or_return
	unmarshal_table(&p, root_table_obj, .EOF, "File") or_return
	return nil
}

// logs the error using log.error, up to the user to call if they want.
log_error :: proc(err_union: Unmarshal_Error) {
	loc_str :: proc(l: Location) -> string {
		@(static) buf: [64]byte
		return fmt.bprintf(buf[:], "%s:%d:%d", l.path, l.row, l.col)
	}


	switch err in err_union {
	case Unexpected_Char:
		if exp, did_expect := err.expected.?; did_expect {
			log.errorf("%s Unexpected character %c expected %c", loc_str(err.location), err.got, exp)
		} else {
			log.errorf("%s Unexpected character %c", loc_str(err.location), err.got)
		}
	case Unexpected_Token:
		if exp, did_expect := err.expected.?; did_expect {
			log.errorf("%s Unexpected token %s expected %s", loc_str(err.got.location), token_repr(err.got), exp)
		} else {
			log.errorf("%s Unexpected token %s", loc_str(err.got.location), token_repr(err.got))
		}
	case Type_Error:
		log.errorf(
			"%s Value %s provided did not match expected type %v",
			loc_str(err.got.location),
			token_repr(err.got),
			type_info_of(err.expected), 
		)
	case Missing_Key:
		log.errorf("Schema expected a key %s in table %s", err.key, err.table)
	case Duplicate_Key:
		log.errorf("Key %s was defined twice in table %s", err.key, err.table)
	case Unknown_Key:
		log.errorf("Schema was not expecting key %s in table %s", err.key, err.table)
	case Unsupported_Type:
		// by default I want to panic on this, because it's more likely a bug then a user error
		log.panicf("Unsupported type in schema %v", type_info_of(err.tid))
	case mem.Allocator_Error:
		log.fatal("Unable to allocate memory.")
	}
}

// If 64 "live" items in an expression is less than your use case needs, increase the limit,
// or write a custom evaluator with a dynamic array instead.
Eval_Stack :: sa.Small_Array(64, f64)
//bool indicates if the function was found or not.
Host_Fn :: #type proc(name: string, stack: ^Eval_Stack) -> bool
// default evaluator, one can always build a custom one, up to user to call after unmarshalling.
evaluate :: proc(expr: Expression, host_vars: map[string]f64, host_fn: Host_Fn) -> (f: f64, ok: bool) {
	stack := Eval_Stack{}
	for op_union in expr {
		switch op in op_union {
		case f64:
			sa.push_back(&stack, op)
		case Expr_Op:
			b := sa.pop_back(&stack)
			a := sa.pop_back(&stack)
			switch op {
			case .Add:
				sa.push_back(&stack, a + b)
			case .Multiply:
				sa.push_back(&stack, a * b)
			case .Subtract:
				sa.push_back(&stack, a - b)
			case .Divide:
				sa.push_back(&stack, a / b)
			}
		case Host_Variable:
			if var, exists := host_vars[string(op)]; exists {
				sa.push_back(&stack, var)
			} else {
				log.errorf("Cannot use variable %s, the host has not defined it.", string(op))
				return 0, false
			}
		case Host_Function:
			host_fn(string(op), &stack) or_return
		}
	}
	return sa.pop_back(&stack), true
}

MAX_FIELDS_PER_STRUCT :: 64

//Note to self: Max fields is defined so that we can use a small array container.
//              A full dynamic array can't be deleted when backed by an arena.
//              Likely not a big deal but for some reason it didn't sit right with me.
@(private)
unmarshal_table :: proc(p: ^Parser, v: any, delimiter: Token_Kind, table_name: string) -> Unmarshal_Error {
	v := v
	ti := reflect.type_info_base(type_info_of(v.id))

	if !reflect.is_struct(ti) do return Unsupported_Type{ti.id}
	fields := reflect.struct_fields_zipped(ti.id)

	fields_filled: sa.Small_Array(MAX_FIELDS_PER_STRUCT, i8)
	entry_loop: for p.token.kind != delimiter {
		key := parse_key(p) or_return
		found_idx: i8 = -1
		for field, i in fields {
			if field.name == key {
				if found_idx != -1 do return Duplicate_Key{key, table_name}
				found_idx = i8(i)
			}
		}
		if found_idx == -1 do return Unknown_Key{key, table_name}
		sa.push_back(&fields_filled, found_idx)

		offset := fields[found_idx].offset
		type := fields[found_idx].type

		field_ptr := rawptr(uintptr(v.data) + offset)
		field := any{field_ptr, type.id}
		unmarshal_value(p, field, key) or_return

		//optional comma between key-value pairs.
		if p.token.kind == .Comma do parser_advance(p) or_return
		if p.token.kind == .EOF do break
	}
	_ = parser_expect(p, delimiter) or_return

	for field, i in fields {
		if slice.contains(sa.slice(&fields_filled), i8(i)) do continue
		// At this point the field was not filled, meaning the key was not in the config.
		// if the type is Maybe(T), we can ignore the missing key.
		ti := reflect.type_info_base(field.type)
		if !reflect.is_union(ti) do return Missing_Key{field.name, table_name}
		union_info := ti.variant.(reflect.Type_Info_Union)
		//heres the check for Maybe, could be more robust I think.
		if len(union_info.variants) > 1 do return Missing_Key{field.name, table_name}
	}
	return nil
}

//Note to self: key name is just used for better errors.
@(private)
unmarshal_value :: proc(p: ^Parser, v: any, key_name: string) -> Unmarshal_Error {
	v := v
	tk := p.token
	ti := reflect.type_info_base(type_info_of(v.id))
	type_err := Type_Error{tk, v.id}

	if u, is_union := ti.variant.(reflect.Type_Info_Union); is_union {
		return unmarshal_union(p, v, u, type_err, key_name)
	}

	parser_advance(p) or_return

	#partial switch tk.kind {
	case .Lit_Float:
		if !assign_float(v, tk.data.(f64)) do return type_err
	case .True, .False:
		if !assign_bool(v, tk.kind == .True) do return type_err
	case .Lit_String:
		switch &dst in v {
		case string:
			dst = strings.clone(tk.data.(string)) or_return
		case:
			return type_err
		}
	case .Lit_Int:
		if assign_int(v, tk.data.(int)) do return nil
		if assign_float(v, tk.data.(int)) do return nil
		return type_err
	case .Identifier:
		if !reflect.is_enum(type_info_of(v.id)) do return type_err
		if val, ok := reflect.enum_from_name_any(v.id, tk.data.(string)); ok {
			assert(assign_int(v, val), "Bug: This should not fail given previous checks")
			return nil
		}
		return type_err
	case .Lbracket:
		#partial switch inner_info in ti.variant {
		case reflect.Type_Info_Array:
			unmarshal_array(p, v, inner_info, key_name) or_return
		case reflect.Type_Info_Slice:
			unmarshal_slice(p, v, inner_info, key_name) or_return
		case:
			return type_err
		}
		_ = parser_expect(p, .Rbracket) or_return
	case .Lbrace:
		return unmarshal_table(p, v, .Rbrace, key_name)
	case .Langle_Bracket:
		ops := parse_expression(p) or_return
		if !reflect.is_slice(type_info_of(v.id)) do return type_err
		mem.copy(v.data, &ops, size_of(ops))
	case:
		return Unexpected_Token{got = tk}
	}
	return nil

	// this is ripped almost verbatim from encoding/json.
	unmarshal_union :: proc(
		p: ^Parser,
		v: any,
		u: reflect.Type_Info_Union,
		type_err: Type_Error,
		key_name: string,
	) -> Unmarshal_Error {
		v := v
		ti := reflect.type_info_base(type_info_of(v.id))
		err: Unmarshal_Error
		for variant, i in u.variants {
			variant_any := any{v.data, variant.id}
			variant_p := p^
			if err = unmarshal_value(&variant_p, variant_any, key_name); err != nil do continue

			p^ = variant_p
			raw_tag := i
			if !u.no_nil do raw_tag += 1
			tag := any{rawptr(uintptr(v.data) + u.tag_offset), u.tag_type.id}
			if !assign_int(tag, int(raw_tag)) do return type_err
			return nil
		}
		return err
	}

	unmarshal_array :: proc(
		p: ^Parser,
		v: any,
		arr_info: reflect.Type_Info_Array,
		key_name: string,
	) -> Unmarshal_Error {
		for i in 0 ..< arr_info.count {
			value_any := any{rawptr(uintptr(v.data) + uintptr(i * arr_info.elem_size)), arr_info.elem.id}
			unmarshal_value(p, value_any, key_name) or_return
			//optional final comma
			if i != arr_info.count - 1 {
				_ = parser_expect(p, .Comma) or_return
			} else {
				if p.token.kind == .Comma do parser_advance(p) or_return
			}
		}
		return nil
	}

	unmarshal_slice :: proc(
		p: ^Parser,
		v: any,
		sl_info: reflect.Type_Info_Slice,
		key_name: string,
	) -> Unmarshal_Error {
		backing: []byte
		i: int
		for i = 0; p.token.kind != .Rbracket; i += 1 {
			backing = mem.resize_bytes(backing, (i + 1) * sl_info.elem_size) or_return
			value_any := any{rawptr(uintptr(raw_data(backing)) + uintptr(i * sl_info.elem_size)), sl_info.elem.id}
			unmarshal_value(p, value_any, key_name) or_return
			if p.token.kind == .Comma do parser_advance(p) or_return
		}
		//there is likely a better way to assign a slice but this works for now.
		//Note to self: the size_of([]int) is just size of a slice, doesnt matter slice of what
		mem.copy(v.data, &struct {
				ptr: rawptr,
				len: int,
			}{raw_data(backing), i}, size_of([]int))
		return nil
	}
}

// these are for assigning to struct fields, cant really do x.y = b, when we don't know the type of x or y at compile time.
// Note to self: the $T notation just infers the type for a generic,
//              meaning you don't need to pass the type explicitly in to a generic function.
//              https://odin-lang.org/docs/overview/#implicit-parametric-polymorphism
@(private)
assign_int :: proc(val: any, i: $T) -> bool {
	v := reflect.any_core(val)
	switch &dst in v {
	case i64:
		dst = i64(i)
	case u64:
		dst = u64(i)
	case int:
		dst = int(i)
	case uint:
		dst = uint(i)
	case u8:
		dst = u8(i)
	case i8:
		dst = i8(i)
	case:
		return false
	}
	return true
}

@(private)
assign_float :: proc(val: any, f: $T) -> bool {
	v := reflect.any_core(val)
	switch &dst in v {
	case f32:
		dst = f32(f)
	case f64:
		dst = f64(f)
	case:
		return false
	}
	return true
}

@(private)
assign_bool :: proc(val: any, b: $T) -> bool {
	v := reflect.any_core(val)
	switch &dst in v {
	case bool:
		dst = bool(b)
	case:
		return false
	}
	return true
}
