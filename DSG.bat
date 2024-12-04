@echo off
rem �������
echo ���ڱ��� my_program
g++ -fexec-charset=UTF-8 -o my_program main.cpp graph.cpp DSG.cpp -std=c++17
if %errorlevel% neq 0 (
    echo ����ʧ�ܡ�����Դ�ļ��ͱ��������á�
    pause
    exit /b
) else (
    echo ����ɹ�
)

rem ����ִ���ļ��Ƿ����
if not exist my_program.exe (
    echo �Ҳ��� my_program.exe�������ļ��Ƿ���ڡ�
    pause
    exit /b
) else (
    echo my_program.exe �ļ�����
)

rem ����ִ�� my_program.exe ����

rem ����1
echo ����ִ������1: Facebook_edges.txt
my_program.exe ../dataset/Facebook_edges.txt 20 200
if %errorlevel% neq 0 (
    echo ����1ִ��ʧ��
) else (
    echo ����1ִ�гɹ�
)

rem ����2
echo ����ִ������2: email-Enron_edges.txt
my_program.exe ../dataset/email-Enron_edges.txt 20 200
if %errorlevel% neq 0 (
    echo ����2ִ��ʧ��
) else (
    echo ����2ִ�гɹ�
)

rem ����3
echo ����ִ������3: Brightkite_edges.txt
my_program.exe ../dataset/Brightkite_edges.txt 20 200
if %errorlevel% neq 0 (
    echo ����3ִ��ʧ��
) else (
    echo ����3ִ�гɹ�
)

rem ����4
echo ����ִ������4: Gowalla_edges.txt
my_program.exe ../dataset/Gowalla_edges.txt 20 200
if %errorlevel% neq 0 (
    echo ����4ִ��ʧ��
) else (
    echo ����4ִ�гɹ�
)

rem ����5
echo ����ִ������5: Twitter_edges.txt
my_program.exe ../dataset/Twitter_edges.txt 20 200
if %errorlevel% neq 0 (
    echo ����5ִ��ʧ��
) else (
    echo ����5ִ�гɹ�
)

rem ����6
echo ����ִ������6: Stanford_edges.txt
my_program.exe ../dataset/Stanford_edges.txt 20 200
if %errorlevel% neq 0 (
    echo ����6ִ��ʧ��
) else (
    echo ����6ִ�гɹ�
)

rem ����7
echo ����ִ������7: Google_edges.txt
my_program.exe ../dataset/Google_edges.txt 20 200
if %errorlevel% neq 0 (
    echo ����7ִ��ʧ��
) else (
    echo ����7ִ�гɹ�
)

rem ����8
echo ����ִ������8: Youtube_edges.txt
my_program.exe ../dataset/Youtube_edges.txt 20 200
if %errorlevel% neq 0 (
    echo ����8ִ��ʧ��
) else (
    echo ����8ִ�гɹ�
)

rem ����9
echo ����ִ������9: Baidubaike_edges.txt
my_program.exe ../dataset/Baidubaike_edges.txt 20 200
if %errorlevel% neq 0 (
    echo ����9ִ��ʧ��
) else (
    echo ����9ִ�гɹ�
)

rem ����10
echo ����ִ������10: as-Skitter.txt
my_program.exe ../dataset/as-Skitter.txt 20 200
if %errorlevel% neq 0 (
    echo ����10ִ��ʧ��
) else (
    echo ����10ִ�гɹ�
)

rem ����11
echo ����ִ������11: socfb-konect.txt
my_program.exe ../dataset/socfb-konect.txt 20 200
if %errorlevel% neq 0 (
    echo ����11ִ��ʧ��
) else (
    echo ����11ִ�гɹ�
)

echo ���������
pause
