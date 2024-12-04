@echo off
rem 编译程序
echo 正在编译 my_program
g++ -fexec-charset=UTF-8 -o my_program main.cpp graph.cpp DSG.cpp -std=c++17
if %errorlevel% neq 0 (
    echo 编译失败。请检查源文件和编译器设置。
    pause
    exit /b
) else (
    echo 编译成功
)

rem 检查可执行文件是否存在
if not exist my_program.exe (
    echo 找不到 my_program.exe。请检查文件是否存在。
    pause
    exit /b
) else (
    echo my_program.exe 文件存在
)

rem 批量执行 my_program.exe 命令

rem 命令1
echo 正在执行命令1: Facebook_edges.txt
my_program.exe ../dataset/Facebook_edges.txt 20 200
if %errorlevel% neq 0 (
    echo 命令1执行失败
) else (
    echo 命令1执行成功
)

rem 命令2
echo 正在执行命令2: email-Enron_edges.txt
my_program.exe ../dataset/email-Enron_edges.txt 20 200
if %errorlevel% neq 0 (
    echo 命令2执行失败
) else (
    echo 命令2执行成功
)

rem 命令3
echo 正在执行命令3: Brightkite_edges.txt
my_program.exe ../dataset/Brightkite_edges.txt 20 200
if %errorlevel% neq 0 (
    echo 命令3执行失败
) else (
    echo 命令3执行成功
)

rem 命令4
echo 正在执行命令4: Gowalla_edges.txt
my_program.exe ../dataset/Gowalla_edges.txt 20 200
if %errorlevel% neq 0 (
    echo 命令4执行失败
) else (
    echo 命令4执行成功
)

rem 命令5
echo 正在执行命令5: Twitter_edges.txt
my_program.exe ../dataset/Twitter_edges.txt 20 200
if %errorlevel% neq 0 (
    echo 命令5执行失败
) else (
    echo 命令5执行成功
)

rem 命令6
echo 正在执行命令6: Stanford_edges.txt
my_program.exe ../dataset/Stanford_edges.txt 20 200
if %errorlevel% neq 0 (
    echo 命令6执行失败
) else (
    echo 命令6执行成功
)

rem 命令7
echo 正在执行命令7: Google_edges.txt
my_program.exe ../dataset/Google_edges.txt 20 200
if %errorlevel% neq 0 (
    echo 命令7执行失败
) else (
    echo 命令7执行成功
)

rem 命令8
echo 正在执行命令8: Youtube_edges.txt
my_program.exe ../dataset/Youtube_edges.txt 20 200
if %errorlevel% neq 0 (
    echo 命令8执行失败
) else (
    echo 命令8执行成功
)

rem 命令9
echo 正在执行命令9: Baidubaike_edges.txt
my_program.exe ../dataset/Baidubaike_edges.txt 20 200
if %errorlevel% neq 0 (
    echo 命令9执行失败
) else (
    echo 命令9执行成功
)

rem 命令10
echo 正在执行命令10: as-Skitter.txt
my_program.exe ../dataset/as-Skitter.txt 20 200
if %errorlevel% neq 0 (
    echo 命令10执行失败
) else (
    echo 命令10执行成功
)

rem 命令11
echo 正在执行命令11: socfb-konect.txt
my_program.exe ../dataset/socfb-konect.txt 20 200
if %errorlevel% neq 0 (
    echo 命令11执行失败
) else (
    echo 命令11执行成功
)

echo 批处理完成
pause
