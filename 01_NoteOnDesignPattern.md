# Notes of `Hands-On Design Patterns and Best Practices with Julia`

[TOC]

## 1. Design patterns

- **SOLID**: Single Responsibility, Open/Closed, Liskov Substitution, Interface
- **DRY**:Don't Repeat Yourself
- **KISS**: Keep It Simple, Stupid!
- **POLA**: Principle of Least Astonishment
- **YAGNI**: You Aren't Gonna Need It
-  **POLP**: Principle of Least Privilege

> **Liskov Substitution**: 类的继承
> **Dependency Inversion**: 高层次的模块不应该依赖于低层次的模块，两者都应该依赖于抽象接口；抽象接口不应该依赖于具体实现。而具体实现则应该依赖于抽象接口。

### Reusability

- top-down
- bottom-up

### Performance

- Functions are small and simple
- Numeric data in contiguous memory space (full use of CPU)
- Minimum memory allocation

### Maintainability

- No unused code
- No duplicate code
- Code is concise and short
- Code is clear and simple
- Every function has a single purpose
- Functions within a module work with each other

### Safety

- Minimize set of types, functions and variables
- call function with arguments
- documented return value of a function
- handle missing data
- small scope variables
- handle exceptions

## 2. Julia base

### namespace, modules, packages

#### Create your package

organize code in Julia packages: **PkgTemplates**


```julia
 set the global git user.name and user.email info in the bash before use julia PkgTemplates
 git config --global user.email "songtaogui@sina.com"
 git config --global user.name "songtaogui"

using PkgTemplates

 set pkg templates
template = Template(; license = "MIT", user = "songtaogui", dir="./")

 generate package named TestJL
generate(template, "TestJL")

 will generate a git local repo of "TestJL", with MIT LICENSE:
 ├── LICENSE
 ├── Manifest.toml
 ├── Project.toml
 ├── README.md
 ├── src
 │   └── TestJL.jl
 └── test
     └── runtests.jl

```

#### create your module:

```julia
 creat a module named Calculator
module Calculator

 export functions so that the functions can be loaded when load the module with "using" keyword
export interest, rate

 write functions with comments
"""
    interest(amount, rate)
Calculate interest from an `amount` and interest rate of `rate`.
"""
function interest(amount, rate)
    return amount * (1 + rate)
end

"""
    rate(amount, interest)
Calculate interest rate based on an `amount` and `interest`.
"""
function rate(amount, interest)
    return interest / amount
end

end # module
```

> **Diff: using vs import**
> ![Acrobat_2021-04-27_13-39-00](https://raw.githubusercontent.com/songtaogui/imgbed/master/D%3A%5CWorks%5C77_MyGit%5CimgbedAcrobat_2021-04-27_13-39-00.png)
> **NOTE:**
> - use the using statement when you are using the functionality, but choose the import statement when you need to extend the functionality from the module.
> - bring specific names into the current namespace (只把module中你需要的function进行using)
> - 特殊情况下，如果有重名的function， 就用import，然后用全名调用函数：Calculator.rate；Rater.rate

#### Submodules

用`include`加载submodules
include的顺序就是加载的顺序，所以submodules如果有依赖，顺序不要搞错
```julia
module Calculator

include("SubA.jl")
include("SubB,jl")

end
```

#### 版本控制：semantic versioning scheme

```julia
<major>.<minor>.<patch>-<pre-release>+<build>
```

添加版本依赖信息到package：

- 在package目录打开julia REPL
- "]" 进入Pkg模式
- `activate .` 激活当前目录的project environment
- `add` 添加依赖的包

**Project.toml**: 当前project的唯一标识信息

```julia
name = "Calculator"
uuid = "db3c262e-fe57-55aa-ad54-d3c64a4d0c83"
authors = ["Tom Kwong <tk3369@gmail.com>"]
version = "0.1.0"

[deps]
SaferIntegers = "88634af6-177f-5301-88b8-7819386cfa38"

[extras]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[targets]
test = ["Test"]

```

**Manifest.toml**: 依赖关系

```julia
 This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[CSTParser]]
deps = ["LibGit2", "Test", "Tokenize"]
git-tree-sha1 = "437c93bc191cd55957b3f8dee7794b6131997c56"
uuid = "00ebfdb7-1f24-5e51-bd34-a7502290713f"
version = "0.5.2"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "84aa74986c5b9b898b0d1acaf3258741ee64754f"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "2.1.0"

[[DataStructures]]
deps = ["InteractiveUtils", "OrderedCollections", "Random", "Serialization", "Test"]
git-tree-sha1 = "ca971f03e146cf144a9e2f2ce59674f5bf0e8038"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.15.0"
 ...
 see https://github.com/PacktPublishing/Hands-on-Design-Patterns-and-Best-Practices-with-Julia/blob/master/Chapter02/Calculator/Manifest.toml
```

NOTE: 版本控制的语法，略

### TYPES

**abstract type**:

```julia
abstract type Asset end
 "<:" --> is a subtype of 
abstract type Property <: Asset end

 subtypes() --> list subtypes of a type
 supertype()
```

**composite type**


```julia
 struct keyword
struct Stock <: Equity
    symbol::String
    name::String
end
```

**Mutability**


```julia
mutable struct Stock <: Equity
    symbol::String
    name::String
end
```

NOTE: immuatble objects are usually the better choice：

- mutable 不容易维护，但是会省内存

**Union types, Nothing**


```julia
Union{Int64, BigInt}
Union{T,Nothing}
```

**`isa` vs `<:`**

```julia
isa -> value vs type
<:  -> type  vs type
```

#### Type conversion

`convert()` function

auto conversion:

1. assign value to array
2. assign value to a field of an object
3. Construct an object with `new`
4. assign value to a variable with declared type
5. function with declared return type
6. pass a value to `ccall`

## 3. Functions and Interfaces

### miscellaneous

- 建议函数名用“_”分隔单词
- 数值改动的函数用"!"
- Duck typing: If it walks like a duck, swims like a duck, and quacks like a duck, then it probably is a duck.
- 可选参数, 关键字参数->';'后边, 可变参数->'...'
- splatting
- anonymous function: `x -> x+1`
- `do...end`

    ```julia
    fire(spaceship) do s
        move_up!(s, 100)
        println(s, " launched missile!)
        move_down!(s, 100)
    end
    # do...end block was turned into anonymous function
    ```

### dispatch

**Multiple Dispatch**: julia matched to the narrowest types

**dynamic dispatch**: the proper method is dispatched right at the argument being passed to the meth

#### detect ambiguity dispatch with `test`

```julia
using Test
detect_ambiguities(Main.Foo)
```

### type parameters


```julia
function explode(things::AbstractVector{T}) where {T <: Thing}
    for t in things
        println("Exploding thing => ", t)
    end
end

# enforce type consistency
function group_same_thing(A::T, B::T) where {T <: Thing}
    println("Grouped", A, " and ", B)
end
```

### extract type info

```julia
eltype([1,2.0])
eltype([1, 2.0, "ABC"])
```

### interface

Interfaces are behavioral contracts


```julia
module Vehicle
# 1. Export/Imports
# 2. Interface documentation
# 3. Generic definitions for the interface
# 4. Game logic
end # module
```

### Soft contract


```julia
# soft contracts
engage_wheels!(args...) = nothing
```

### trait

whether a data type exhibits certain behavior

```julia
# trait
has_wheels(vehicle) = error("Not implemented.") # vehicle without wheels
```

## 4. Macros and Metaprogramming

metaprogramming: writing code that generates code

### need for metaprogramming

- more concisely and comprehend
- reduce dev time
- improve performance

e.g.: `@time`:


```julia
# @time equals to:
function timeit(func)
    t1 = now()
    result = func()
    t2 = now()
    elapsed = t2 - t1
    println("It took ", elapsed)
    result
end
```

`@unroll`: unroll loops:


```julia
using Unrolled

@unroll function hello(xs)
    @unroll for i in xs
        println("hello: ", i)
    end
end

function hello_loop(xs)
    for i in xs
        println("hello_loop: ", i)
    end
end

# unroll the for loop:
seq=tuple(1:3...)

@code_lowered(hello(seq))
# CodeInfo(
# 1 ─│ %1  = Base.getindex(xs, 1)
# │  │       i@_3 = %1
# │  │       Main.println("hello: ", i@_3)
# │  │ %4  = Base.getindex(xs, 2)
# │  │       i@_4 = %4
# │  │       Main.println("hello: ", i@_4)
# │  │ %7  = Base.getindex(xs, 3)
# │  │       i@_5 = %7
# │  │       Main.println("hello: ", i@_5)
# │  │ %10 = Main.nothing
# └──│       return %10
#    └
# )

@code_lowered(hello_loop(seq))
# CodeInfo(
# 1 ─ %1  = xs
# │         @_3 = Base.iterate(%1)
# │   %3  = @_3 === nothing
# │   %4  = Base.not_int(%3)
# └──       goto #4 if not %4
# 2 ┄ %6  = @_3
# │         i = Core.getfield(%6, 1)
# │   %8  = Core.getfield(%6, 2)
# │         Main.println("hello_loop: ", i)
# │         @_3 = Base.iterate(%1, %8)
# │   %11 = @_3 === nothing
# │   %12 = Base.not_int(%11)
# └──       goto #4 if not %12
# 3 ─       goto #2
# 4 ┄       return nothing
# )
```

### expressions is just a data structure in julia

- abstract syntax tree (AST)

    `Meta.parse`: parse a string into an _Expr_ object
    `dump`: print the structure of an Expr


    ```julia
    Meta.parse("x + y") |> dump
    # Expr
    #   head: Symbol call
    #   args: Array{Any}((3,))
    #     1: Symbol +
    #     2: Symbol x
    #     3: Symbol y
    ```
    ![Acrobat_2021-04-28_20-12-10](https://raw.githubusercontent.com/songtaogui/imgbed/master/D%3A%5CWorks%5C77_MyGit%5CimgbedAcrobat_2021-04-28_20-12-10.png)

- manually expr

    ```julia
    Expr(:call, :sin, Expr(:call, :+, :x, :y))
    # :(sin(x + y))

    ex = :(sin(x + y))
    # :(sin(x + y))
    ```

    **NOTE:**:
    1. `:( expr )` + `begin ... end` block: work with multiple expressions:

        ```julia
        :(begin
        x = 1
        y = 2
        end)
        ```
    2. `quote ... end` can also create expr:

        ```julia
        quote
            x = 1
            y = 2
        end
        ```

    3. `eval()`: evaluation the expression **in the current module**

    4. _interpolation_ (插值) and _splatting_ : dynamically create expressions

        ```julia
        julia> v = [1,2,3]
        julia> quote
                    max($(v)...) #插值的展开用()控制
            end
        # quote
        #     #= REPL[39]:2 =#
        #     max([1, 2, 3]...)
        # end

        julia> quote
                    max($(v...)) #插值的展开用()控制
            end
        # quote
        #     #= REPL[40]:2 =#
        #     max(1, 2, 3)
        # end
        ```
    5. `QuoteNode()`: assign an actual symbol to a variable

#### Macro

**Why Macro? Macro vs functions:**

> - Macro expansion happens during compilation
> - expression from a macro can be executed within the current scope rather than the global scope

##### write Macro with `macro`

NOTE: macro must return expressions


```julia
macro hello()
    return :(
    for i in 1:3
        println("hello world")
    end
    )
end

@hello()
# macros can be called without parentheses:
@hello
```

##### pass arguments to macro


```julia
# pass literals:
macro hello(n)
    return :(
    for i in 1:$n
        println("hello world")
    end
    )
end

@hello(2) # or: @hello 2

# macro arguments pass expressions
# function arguments pass values
```

##### macro expansion

`@macroexpand`: see the expression of a macro without running it

- any macro is expanded when invoked in the REPL
- macro is expanded as part of the function definition process

    > NOTE: it is important to redefine the function again if any of the macros being used have been changed.

- use `esc` to place the expression directly in the syntax tree without letting the compiler resolve it

    ```julia
    # a wrong example:
    # --------------------------
    macro squared1(ex)
        return :($(ex) * $(ex))
    end

    function foo1()
        x=2
        return @squared1 x
    end
    # --------------------------
    julia> foo()
    ERROR: UndefVarError: x not defined
    # WHY?
    # @square 在函数定义时直接展开了，用的是当前
    # module（main）中的x，而函数中定义的是local x:
    julia> @code_lowered foo1()
    CodeInfo(
    1 ─      x = 2
    │   %2 = Main.x * Main.x
    └──      return %2
    )

    # 用 esc() 强制把表达式直接放到语法树中：
    # Prevents the macro hygiene pass from
    # turning embedded variables into gensym
    # variables

    # a better example:
    # --------------------------
    macro squared2(ex)
        return :($(esc(ex)) * $(esc(ex)))
    end

    function foo2()
        x=2
        return @squared2 x
    end
    # --------------------------
    julia> @code_lowered foo2()
    CodeInfo(
        1 ─      x = 2
    │   %2 = x * x
    └──      return %2
    # --------------------------
    julia> foo2()
    4
    ```

### developing macros

- manipulating expressions through AST

    ```julia
    # design a macro that takes a simple funcction call expression and calls the same function again with the result:
    # @compose_twice sin(x) -> sin(sin(x))

    macro compose_twice(ex)
        @assert ex.head == :call
        @assert length(ex.args) == 2
        me = copy(ex)
        ex.args[2] = me
        return ex
    end

    # @assert: 类似perl中的die；这里用来确保输入的function只有一个参数；
    # 嵌套一次的话， 就直接把AST的参数替换成原来的AST就行了
    ```

- macro hygiene

    避免宏与其他内容冲突：


    ```julia
    # @ntimes -> 重复一个函数n次
    macro ntimes(n, ex)
        quote
            times = $(esc(n))
            for i in 1:times
                $(esc(ex))
            end
        end
    end

    # 如果 times 同时在宏外被赋值，是否会冲突？
    function foo()
        times = 0
        @ntimes 3 println("hello world")
        println("times = ", times)
    end

    # 不冲突！
    julia> foo()
    hello world
    hello world
    hello world
    times = 0

    # WHY？
    # Julia的宏在展开时把变量都进行了重命名，确保不会有冲突：

    julia> @macroexpand(@ntimes 3 println("hello world"))
    quote
        #= REPL[98]:4 =#
        var"#58#times" = 3  # times 被重命名成 #58#times
        #= REPL[98]:5 =#
        for var"#59#i" = 1:var"#58#times"  # i 变成了 #59#i
            #= REPL[98]:6 =#
            println("hello world")
        end
    end

    # nice!
    ```

- nonstandard string literals

    定义`非标准字符串语义`, julia的regex就是一种:

    ```julia
    macro r_str(p)
        Regex(p)
    end
    # XXX_str 定义宏， 然后就可以用 XXX"non_std string"的格式去调用这个宏了
    #
    # 用 r"regex_patterns" 而不是Regex("regex_patterns") 的好处是：
    #
    # r"regex_patterns" 调用宏，所以第一次出现时就直接展开到AST中了，后边不用重复编译，可以直接放到循环中。而用Regex()的话，因为是函数，如果直接放到循环中，就每次都要重新编译，所以直接放到循环中就会降低性能，要在循环外赋值给变量，然后在循环内调用变量。

    for line = lines
        m = match(r"^abc$", line)
        if m === nothing
            # nothing
        else
            # something
        end
    end

    # 用 Regex() 就要：
    re = Regex("^abc$")
    for line = lines
        m = match(re, line)
        if m === nothing
            # nothing
        else
            # sth
        end
    end

    # 思考：
    # 如果循环内的 Regex 每次都变化的话, 那把Regex放到循环内跟宏应该差不多效果。
    # Regex如果依赖变量变化，用Regex()函数形式更容易实现，可以直接插值，'$'需要转义。
    # 总之，如果在循环内用的话，直接用非标准字符串的形式更好。
    ```

### using generated functions

**What Why and How of generated functions**

> - what: 生成函数能根据其参数类型生成专用代码
> - why: 解决Macro无法考虑到类型的问题，更专业地：[一种兼顾简便, 编程效率和性能的staging实现](https://thautwarm.github.io/Site-32/Design/Staging%E5%92%8CJulia%E7%94%9F%E6%88%90%E5%87%BD%E6%95%B0.html)
> - how: `@generate`

**Macro, generated function, and function:**

> Macro直接编辑AST，在传递类型之前就编译了，参数是表达式， 返回表达式
> generated function 的参数是类型，返回表达式
> function的参数是变量的值，返回运行结果


```julia
@generated function doubled(x)
    @show x
    return :(2 * x)
end

# 第一次运行doubled(2)的时候，根据2的类型，生成doubled()函数
julia> doubled(2)
x = Int64
4
# 后边再doubled(2)就是直接运行已经编译的同类型函数了
julia> doubled(3)
6

# 再运行一个新类型的时候会重新编译一个新的doubled()函数：
julia> doubled(3.0)
x = Float64
6.0

julia> doubled(4.0)
8.0

# 假设以后有一种新function可以高效计算浮点数：叫double_super_duper
# 用生成函数可以方便地按照类型指定调用的function：
@generated function doubled(x)
    if x <: AbstractFloat
        return :( double_super_duper(x))
    else
        return :(2 * x)
    end
end

```

> **QUST: 多重派发和生成函数的区别？**
> 我的理解：生成函数更顶层，可以生成包含多重派发的function


[拓展阅读：Staging和Julia生成函数](https://thautwarm.github.io/Site-32/Design/Staging%E5%92%8CJulia%E7%94%9F%E6%88%90%E5%87%BD%E6%95%B0.html#stagingjulia)

> - **Staging**: 通过将程序运行时分割成更多的编译/运行期，每一阶段，都会运用更丰富的信息生成更高效的程序

## 5. Reusability Patterns

### 5.1 delegation pattern

delegation: 委托模式

forward: 转发


```julia
# 定义Account数据结构，包含数据存取器accessors和一些方法：deposit withdraw transfer
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Object Definition
mutable struct Account
    account_number::String
    balance::Float64
    date_opened::Date
end

# Accessors
account_number(a::Account) = a.account_number
balance(a::Account) = a.balance
date_opened(a::Account) = a.date_opened

# Functions

function deposit!(a::Account, amount::Real)
    a.balance += amount
    return a.balance
end

function withdraw!(a::Account, amount::Real)
    a.balance -= amount
    return a.balance
end

function transfer!(from::Account, to::Account, amount::Real)
    # println("Transferring ", amount, " from account ",
    #     account_number(from), " to account ", account_number(to))
    withdraw!(from, amount)
    deposit!(to, amount)
    return amount
end
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 定义一个新的结构，转发Account中的accessor和方法
"""
`SavingsAccount` is a kind of `Account` that generates interest daily.
We will use Delegation pattern to gain existing `Account` functionality.
"""
struct SavingsAccount
    acct::Account
    interest_rate::Float64
    # 创建新的Account结构，加上interest_rate
    SavingsAccount(account_number, balance, date_opened, interest_rate) = new(
        Account(account_number, balance, date_opened),
        interest_rate
    )
end

# Forward accessors 转发存取器
account_number(sa::SavingsAccount) = account_number(sa.acct)
balance(sa::SavingsAccount) = balance(sa.acct)
date_opened(sa::SavingsAccount) = date_opened(sa.acct)

deposit!(sa::SavingsAccount, amount::Real) = deposit!(sa.acct, amount)
withdraw!(sa::SavingsAccount, amount::Real) = withdraw!(sa.acct, amount)
transfer!(sa1::SavingsAccount, sa2::SavingsAccount, amount::Real) = transfer!(
    sa1.acct, sa2.acct, amount)

# new accessor
interest_rate(sa::SavingsAccount) = sa.interest_rate

# new behavior
function accrue_daily_interest!(sa::SavingsAccount) 
    interest = balance(sa.acct) * interest_rate(sa) / 365
    deposit!(sa.acct, interest)
end


# 以上转发操作可以通过宏更简单的实现：

using Lazy: @forward

@forward SavingsAccount.acct account_number, balance, date_opened
@forward SavingsAccount.acct deposit!, withdraw!

transfer!(from::SavingsAccount, to::SavingsAccount, amount::Real) = transfer!(
    from.acct, to.acct, amount)

# @forward 对象 [待转发的方法元组]
# @forward 不能转发transfer!函数，因为它有两个参数，所以要手动转发
# 有其他包支持多参数的转发，如TypedDelegation.jl
# QUST:转发和类继承的区别之一：
# 类继承不能DIY哪些属性被继承，转发可以DIY
```







## 6. Performance Patterns

### 6.1 The global constant pattern

- 全局变量方便，但是很慢， 一个调用全局变量的函数在编译的时候，必须要考虑到全局变量会改变，所以要把函数编译成能适用各种类型的模式。
- 用全局常量可以有效避免：

    ```julia
    # global variable is slow!
    variable = 10

    function add_using_global_variable(x)
    return x + variable
    end

    julia> @btime add_using_global_variable(10)
    22.992 ns (0 allocations: 0 bytes)
    20

    # global constant is as fast as local!
    const constant = 10

    function add_using_global_constant(x)
    return constant + x
    end

    julia> @btime add_using_global_constant(10)
    1.300 ns (0 allocations: 0 bytes)
    20
    ```
- 为啥常量就快呢？

    > - 数值不变
    > - 类型不变

    Julia的编译器优化涉及：常量折叠->常量传播->死码删除，最终呈现出的是非常简洁的实现，所以快：

    ```julia
    function constant_folding_example()
        a = 2 * 3
        b = a + 1
        return b > 1 ? 10 : 20
    end
    # 这个函数一直返回10，因为a、b都是常量，判断实际上也没有意义，Julia的解释器懂么？
    julia> @code_typed constant_folding_example()
    CodeInfo(
    1 ─     return 10
    ) => Int64
    # @code_typed宏展示类型解析后的代码逻辑，好家伙，直接就返回10了，Julia好会啊！
    ```

- 如果必须要用到全局变量怎么办？
    -> 在调用的函数内部添加类型注释：

    ```julia
    variable = 10

    function add_using_global_variable_typed(x)
        return x + variable::Int
    end

    julia> @btime add_using_global_variable_typed(10)
    2.700 ns (0 allocations: 0 bytes)
    20

    # 比直接用全局变量快不少，但还是比不上用全局常量
    ```

- 还有更快的方法么？
    -> 有！ 把全局变量当作参数传递给函数：

```julia
function add_by_passing_global_variable(x, v)
    return x + v
end

julia> @btime add_by_passing_global_variable(10, $variable)
  1.400 ns (0 allocations: 0 bytes)
20

# Nice! 跟用常量一样快了！
```

- 还有其他方法？
    -> **global variable place holder**
    只要在编译的时候所有变量的类型确定，Julia就能生成最优化的代码，所以可以创建一个常量占位符，然后把变量存进去：

    ```julia
    # julia中的Ref对象是类型确定的占位符：
    julia> Ref(10)
    Base.RefValue{Int64}(10)

    julia> Ref("abc")
    Base.RefValue{String}("abc")

    # 可以用空索引操作来获取占位符中的值，也可以赋值
    julia> a = Ref("abc")
    Base.RefValue{String}("abc")

    julia> a[]
    "abc"

    julia> a[] = "efg"
    "efg"

    # 用全局占位符看看效率如何：
    const semi_constant = Ref(10)

    function add_using_global_semi_constant(x)
        return x + semi_constant[]
    end

    julia> @btime add_using_global_semi_constant(10)
    3.000 ns (0 allocations: 0 bytes)
    20

    # 也还不错！
    ```