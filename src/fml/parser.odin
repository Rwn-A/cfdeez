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
package fml

import "core:mem"
import "core:strconv"
import "core:strings"
import "core:fmt"

Location :: struct {
	row:  int,
	col:  int,
	path: string,
}


Token_Data :: union {
	int,
	f64,
	string,
}

Token_Kind :: enum {
	Lit_String,
	Lit_Int,
	Lit_Float,
	True,
	False,
	Identifier,
	Colon,
	Comma,
	Lbrace,
	Rbrace,
	Lbracket,
	Rbracket,
	Langle_Bracket,
	Rangle_Bracket,

	//tokens only allowed in expression context
	Expr_Dash,
	Expr_Plus,
	Expr_Lparen,
	Expr_Rparen,
	Expr_Asterisk,
	Expr_SlashForward,
	Expr_Dot,
	EOF,
}


Token :: struct {
	kind:     Token_Kind,
	data:     Token_Data,
	location: Location,
}

// once called, you must use the string before calling again or it will be overwritten
token_repr :: proc(tk: Token, allocator := context.allocator) -> string {
	@(static) buf: [64]byte
	switch data in tk.data{
		case string: return data
		case int: return strconv.itoa(buf[:], data)
		case f64: return strconv.ftoa(buf[:], data, 'f', 8, 64)
		case: return fmt.enum_value_to_string(tk.kind) or_else ""
	}
}

@(private)
Lexer :: struct {
	file_contents: []byte,
	position:      int,
	location:      Location,
	in_expression: bool,
	// I'm not stoked about the lexer being aware of expression context
	// Leaving it to parse time added some difficulty dealing with negative number definitions.
	// For example, x: -123, has to be lexed as a literal int because we can't do expressions in normal fields.
	// But x: <10 - 5>, this would break if both were lexed as literal int.
	// the solution I came up with is to lex the - differently if we are in an expression or not.
}

@(private)
expr_only_chars := strings.ascii_set_make("+-/*()$.") or_else panic("Tried to use non-ascii as expression char.")

@(private)
lexer_next :: proc(l: ^Lexer) -> (tk: Token, err: Unmarshal_Error) {
	EOF :: 0

	for c := lexer_char(l); strings.is_ascii_space(rune(c)); c = lexer_char(l) {
		lexer_advance(l)
		if c != '\n' do continue
		l.location.row += 1;l.location.col = 1
	}

	current_token := Token {
		kind     = .EOF,
		location = l.location,
	}
	c := lexer_char(l)

	if c == EOF do return current_token, nil

	switch c {
	case '{':
		current_token.kind = .Lbrace
	case '}':
		current_token.kind = .Rbrace
	case '[':
		current_token.kind = .Lbracket
	case ']':
		current_token.kind = .Rbracket
	case ':':
		current_token.kind = .Colon
	case ',':
		current_token.kind = .Comma
	case 'a' ..= 'z', 'A' ..= 'Z':
		current_token.kind, current_token.data = lex_ident_or_keyword(l)
	case '-':
		if l.in_expression {
			current_token.kind = .Expr_Dash
		} else {
			current_token.kind, current_token.data = lex_int_or_float(l) or_return
		}
	case '0' ..= '9':
		current_token.kind, current_token.data = lex_int_or_float(l) or_return
	case '"':
		{
			start := l.position
			for p := lexer_peek(l); p != '"'; p = lexer_peek(l) {
				if p == EOF || p == '\n' {
					return {}, Unexpected_Char{p, '"', l.location}
				}
				lexer_advance(l)
			}
			lexer_advance(l) // skip end "
			text := cast(string)l.file_contents[start + 1:l.position]
			current_token.kind = .Lit_String;current_token.data = text
		}
	case '#':
		{
			for lexer_peek(l) != '\n' && lexer_peek(l) != EOF do lexer_advance(l)
			lexer_advance(l)
			return lexer_next(l)
		}
	case '+':
		current_token.kind = .Expr_Plus
	case '*':
		current_token.kind = .Expr_Asterisk
	case '/':
		current_token.kind = .Expr_SlashForward
	case '(':
		current_token.kind = .Expr_Lparen
	case ')':
		current_token.kind = .Expr_Rparen
	case '.':
		current_token.kind = .Expr_Dot
	case '<':
		current_token.kind = .Langle_Bracket
		l.in_expression = true
	case '>':
		current_token.kind = .Rangle_Bracket
		if !l.in_expression {
			return {}, Unexpected_Char{got = c, location = l.location}
		}
		l.in_expression = false
	case:
		return {}, Unexpected_Char{got = c, location = l.location}
	}

	lexer_advance(l)

	if !l.in_expression && strings.ascii_set_contains(expr_only_chars, c) {
		return {}, Unexpected_Char{got = c, location = l.location}
	}

	return current_token, nil

	lex_ident_or_keyword :: proc(l: ^Lexer) -> (Token_Kind, Token_Data) {
		start := l.position
		loop: for {
			switch lexer_peek(l) {
			case 'a' ..= 'z', 'A' ..= 'Z', '0' ..= '9', '_':
				lexer_advance(l)
			case:
				break loop
			}
		}
		text := cast(string)l.file_contents[start:l.position + 1]

		switch text {
		case "true":
			return .True, nil
		case "false":
			return .False, nil
		}

		return .Identifier, text
	}

	lex_int_or_float :: proc(l: ^Lexer) -> (Token_Kind, Token_Data, Unmarshal_Error) {
		start := l.position
		// the - symbol is only allowed after an e, or E, for something like 1e-10.
		allow_dash: bool
		loop: for {
			p := lexer_peek(l)
			switch p {
			case '0' ..= '9', 'x', 'X', '.', '_', 'a' ..= 'f', 'A' ..= 'F':
				allow_dash = false
				if p == 'e' || p == 'E' do allow_dash = true
				lexer_advance(l)
			case:
				if p == '-' && allow_dash {
					lexer_advance(l)
					continue loop
				}
				break loop
			}
		}

		text := cast(string)l.file_contents[start:l.position + 1]
		if val, ok := strconv.parse_int(text); ok {
			return .Lit_Int, val, nil
		}
		if val, ok := strconv.parse_f64(text); ok {
			return .Lit_Float, val, nil
		}


		return nil, nil, Type_Error{Token{.Lit_Float, text, l.location}, f64}
	}

	lexer_char :: #force_inline proc(l: ^Lexer) -> byte {
		return l.file_contents[l.position] if l.position < len(l.file_contents) else 0
	}

	lexer_peek :: #force_inline proc(l: ^Lexer) -> byte {
		return l.file_contents[l.position + 1] if l.position + 1 < len(l.file_contents) else 0
	}

	lexer_advance :: #force_inline proc(l: ^Lexer) {
		l.position += 1;l.location.col += 1
	}
}

@(private)
Parser :: struct {
	lexer: Lexer,
	token: Token,
	peek:  Token,
}

//I often forget to handle the bool here hence the require_results.
@(require_results, private)
parser_advance :: proc(p: ^Parser) -> Unmarshal_Error {
	p.token = p.peek
	p.peek = lexer_next(&p.lexer) or_return
	return nil
}

@(private)
parser_expect :: proc(p: ^Parser, kind: Token_Kind) -> (tk: Token, err: Unmarshal_Error) {
	if p.token.kind != kind {
		return {}, Unexpected_Token{got = p.token, expected = kind}
	}
	tk = p.token
	parser_advance(p) or_return
	return tk, nil
}

@(private)
parser_init :: proc(contents: []byte, filepath: string) -> (p: Parser, err: Unmarshal_Error) {
	p = Parser {
		lexer = Lexer{contents, 0, {1, 1, filepath}, false},
	}
	p.token = lexer_next(&p.lexer) or_return
	p.peek = lexer_next(&p.lexer) or_return

	return p, nil
}

@(private)
parse_key :: proc(p: ^Parser) -> (key: string, err: Unmarshal_Error) {
	tk := parser_expect(p, .Identifier) or_return
	_ = parser_expect(p, .Colon) or_return
	return tk.data.(string), nil
}

Expression :: []Expression_Node
Expression_Node :: union {
	f64,
	Host_Variable,
	Host_Function,
	Expr_Op,
}
Host_Variable :: distinct string
Host_Function :: distinct string
Expr_Op :: enum {
	Add,
	Subtract,
	Multiply,
	Divide,
}

@(private)
Precedence :: enum {
	Lowest,
	Sum,
	Product,
	Unary,
	Highest,
}

@(private)
Token_Kind_Set :: bit_set[Token_Kind]

@(private)
op_tokens: Token_Kind_Set : {.Expr_Plus, .Expr_Dash, .Expr_Asterisk, .Expr_SlashForward}

@(private)
operator_precedence :: proc(kind: Token_Kind) -> Precedence {
	#partial switch kind {
	case .Expr_Plus, .Expr_Dash:
		return .Sum
	case .Expr_Asterisk, .Expr_SlashForward:
		return .Product
	case:
		return .Lowest
	}
}

@(private)
parse_expression :: proc(p: ^Parser) -> (v: Expression, err: Unmarshal_Error) {
	out := make([dynamic]Expression_Node)

	parse_primary_expression(p, .Lowest, &out) or_return
	if p.peek.kind != .Rangle_Bracket {
		return {}, Unexpected_Token{got = p.peek, expected = .Rangle_Bracket}
	}
	parser_advance(p) or_return
	parser_advance(p) or_return
	return out[:], nil

	parse_primary_expression :: proc(
		p: ^Parser,
		prec: Precedence,
		out: ^[dynamic]Expression_Node,
	) -> (
		err: Unmarshal_Error,
	) {
		#partial switch p.token.kind {
		case .Lit_Int:
			append(out, f64(p.token.data.(int)))
		case .Lit_Float:
			append(out, p.token.data.(f64))
		case .Expr_Lparen:
			parse_grouped(p, out) or_return
		case .Identifier:
			if p.peek.kind == .Expr_Lparen {
				parse_host_function(p, out) or_return
			} else {
				name := strings.clone(p.token.data.(string)) or_return
				append(out, Host_Variable(name))
			}
		case .Expr_Dash:
			parse_unary_expr(p, out) or_return
		case:
			return Unexpected_Token{got = p.peek}
		}
		for p.peek.kind != .Rangle_Bracket && prec < operator_precedence(p.peek.kind) {
			if p.peek.kind not_in op_tokens do return nil
			parser_advance(p) or_return
			parse_binary_expr(p, out) or_return
		}
		return
	}

	parse_binary_expr :: proc(p: ^Parser, out: ^[dynamic]Expression_Node) -> (err: Unmarshal_Error) {
		op_token := p.token.kind
		prec := operator_precedence(op_token)

		parser_advance(p) or_return
		parse_primary_expression(p, prec, out) or_return

		#partial switch op_token {
		case .Expr_Plus:
			append(out, Expr_Op.Add)
		case .Expr_Dash:
			append(out, Expr_Op.Subtract)
		case .Expr_Asterisk:
			append(out, Expr_Op.Multiply)
		case .Expr_SlashForward:
			append(out, Expr_Op.Divide)
		case:
			return Unexpected_Token{got = p.peek}
		}
		return
	}

	parse_unary_expr :: proc(p: ^Parser, out: ^[dynamic]Expression_Node) -> (err: Unmarshal_Error) {
		_ = parser_expect(p, .Expr_Dash) or_return
		parse_primary_expression(p, .Lowest, out) or_return
		append(out, -1)
		append(out, Expr_Op.Multiply)
		return
	}

	parse_grouped :: proc(p: ^Parser, out: ^[dynamic]Expression_Node) -> (err: Unmarshal_Error) {
		parser_advance(p) or_return
		parse_primary_expression(p, .Lowest, out) or_return
		if p.peek.kind != .Expr_Rparen {
			return Unexpected_Token{got = p.peek, expected = .Expr_Lparen}
		}
		parser_advance(p) or_return
		return

	}

	parse_host_function :: proc(p: ^Parser, out: ^[dynamic]Expression_Node) -> (err: Unmarshal_Error) {
		name_tk := parser_expect(p, .Identifier) or_return
		parser_advance(p) or_return

		arg_builder := make([dynamic]Expression_Node)
		defer delete(arg_builder)
		for p.token.kind != .Expr_Rparen {
			parse_primary_expression(p, .Lowest, &arg_builder) or_return
			parser_advance(p) or_return
			if p.token.kind == .Comma {
				parser_advance(p) or_return
				continue
			}
		}

		#reverse for arg in arg_builder do append(out, arg)
		name := strings.clone(name_tk.data.(string)) or_return
		append(out, Host_Function(name))

		return nil
	}
}
