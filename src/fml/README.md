# FML
Fluid Markup Language. A configuration file format.

**This README is specifically for the configuration language itself not how it interfaces with the fluid sim.**

FML is similar to JSON. It’s a file of key-value pairs. But unlike JSON:
- Keys are not quoted
- Commas between keys are optional
- Trailing commas between keys and array values are allowed
- Comments are allowed

Also, FML is designed as a **configuration format**, not a data exchange format.

## Why does this exist?
There are some modified JSON specs out there that make JSON a bit friendlier to read and I could add expressions on top.
The main issue is that Odin’s built-in JSON library isn’t nearly as strict as I wanted. So I would have had to roll my own.
Instead of being limited by the JSON spec, I'd rather just make a new format, especially so expressions can have custom syntax.

I also considered TOML, but the third-party TOML parser for Odin doesn’t support unmarshalling into structs.
I’d have to build that myself too.

YAML was never even on the table.

## Syntax
Every thing in FML is a key value pair. `some_key: "some_value"`. Below are all the different types of values you can use.
```yaml
string_key: "string_value"     # some arbitrary text
bool_key: true                 # or false
integer_key: 123               # hexadecimal and scientific notation supported.
float_key: 12.3                # hexadecimal and scientific notation supported.
enum_key: North                # An enum is one variant of a group of specific text.
array_key: [123, 231, 312]     # an array is a group of values, all values must be the same type.
table_key: { inner: "inner" }  # A table is another section of any number of key value pairs
expr_key: <x * 4 + cos(pi) >   # See below.
```

### Expression
The expression is used to define some math the reader of the config should perform. The writer of the config and the reader
should agree on what variables, `x` and `pi` in the example above, are provided. They should also agree on what functions
are used, `cos` in the above example. The idea for these expressions was mainly for defining custom initial conditions and
inflow profiles for my fluid simulation. Perhaps they could be useful outside of that.

## Usage
FML is simple to use if you are familiar with Odin, the below code snippet demonstrates a basic example.
```odin
import fml

import "core:fmt"

main :: proc() {
    my_schema := struct {
        x: int,
        y: string,
        some_expression: fml.Expression,
    }{}

    my_file := "x: 10, y: \"Hello, World\", some_expression: <1 + 3 * cos(r)> "

    err := fml.unmarshal(transmute([]byte)my_file, "some_dummy_filepath", &my_schema, context.temp_allocator)
    defer free_all(context.temp_allocator)

    if err != nil {
        fml.log_error(err)
        return
    }

    fmt.printf("%d, %s, ", my_schema.x, my_schema.y)

    host_fn :: proc(name: string, stack: ^fml.Eval_Stack) -> bool {
        switch name {
        case "cos": sa.push_back(stack, math.cos(sa.pop_back(stack)))
        case: log.errorf("function %s is unknown.", name); return false
        }
        return true
    }

    variables := make(map[string]f64)
    variables["r"] = 2 * 3.14159265

    result, ok := fml.evaluate(schema.some_expression, variables, host_fn)

    fmt.println(result)
}
```
We expect the below code to output `10, Hello, World, 4`.

### Defining a schema
Only a subset of Odin types are currently allowed in FML schema.
```
u8, i8, uint, int
bool
string
f64, f32
arrays, slices, structs and unions
```

There are also two special types.
1. usage of `Maybe(T)` indicates to the unmarshalling code that the field can be ommitted in the config file.
2. usage of fml.Expression indicates to the unmarshalling code that an expression is expected.

