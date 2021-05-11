# notes of Julia basic functions

## 1. file | system | IO ...

Filesystem:

```julia
pwd()
cd()
readdir()
walkdir() # 文件迭代器
mkdir() # mkdir in shell
mkpath() # mkdir -p in shell
symlink()
readlink()
chmod()
chown()
RawFD # Type
stat() #返回数据结构，包含一堆文件信息：

a = stat("note.jl")
a.
blksize blocks   ctime    device   gid
inode    mode     mtime    nlink    rdev
size     uid

lstat() # stat for links
ctime() mtime() # stat(file).ctime stat(file).mtime
filemode() # stat(file).mode
filesize() # stat(file).size
uperm() # get the user permissions of the file
# 01	Execute Permission
# 02	Write Permission
# 04	Read Permission
gperm() # group permissions
operm() # owner permissions
cp()
download()
mv()
rm()
touch()
tempname()
tempdir()
mktemp()
mktempdir()

# 判断
isblockdev() # if path is a block device
ischardev() # if path is a character device
isdir() # if path is a directory
isfile()
ispath()
isfifo()
ismount()
issetgid() # setgid flag:
issetuid() # setuid
issocket() # socket
issticky() # sticky bit set
isabspath()
isdirpath()

homedir()  # ~
dirname() basename()
joinpath()
abspath()
normpath() # rfmt path by removing "." and ".."
realpath() # expanding links
relpath() # relative path
expanduser()
splitdir() # ( "dirname", "filename")
splitdrive() #win
splitext() # get extension
# julia> splitext("abcd/aa.b.c")
# ("abcd/aa.b", ".c")

```

I/O:


```julia
stdin stdout stderr # constants
open()
IOStream # Type
IOBuffer # Type
take!(b:IOBuffer) # 从IOBuffer中取出，buffer之后就回到初始状态

# IOStream
fdio() # 创建IOStream
flush()
close()
write()
read()
read!()
readbytes!()
unsafe_read()
unsafe_write()
readeach()
peek() # 不改变stream中的position
position()
seek()
seekstart() seekend()
skip()
mark() unmark() reset() ismarked()
eof()
isreadonly()
iswritable()
isreadable()
isopen()
fd()
redirect_stdout()
redirect_stderr()
redirect_stdin()
readchomp() # => chomp(read(x, String)) 全加载，只chomp末尾，不是逐行加载
truncate()
skipchars()
countlines()
# julia> countlines(IOBuffer("a\nb\nc"))
# 3

# julia> countlines(IOBuffer("a\nb\nc"),eol = '.')
# 1
PipeBuffer()
readavailable()
IOContext()

# TEST IO
show()
summary()
print()
println() # add new line
printstyled()
```
