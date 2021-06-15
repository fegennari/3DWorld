!include x64.nsh
Name "libjpeg-turbo SDK for Visual C++ 64-bit"
OutFile "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\${BUILDDIR}libjpeg-turbo-2.1.1-vc64.exe"
InstallDir "c:\libjpeg-turbo64"

SetCompressor bzip2

Page directory
Page instfiles

UninstPage uninstConfirm
UninstPage instfiles

Section "libjpeg-turbo SDK for Visual C++ 64-bit (required)"
!ifdef WIN64
	${If} ${RunningX64}
	${DisableX64FSRedirection}
	${Endif}
!endif
	SectionIn RO
!ifdef GCC
	IfFileExists $SYSDIR/libturbojpeg.dll exists 0
!else
	IfFileExists $SYSDIR/turbojpeg.dll exists 0
!endif
	goto notexists
	exists:
!ifdef GCC
	MessageBox MB_OK "An existing version of the libjpeg-turbo SDK for Visual C++ 64-bit is already installed.  Please uninstall it first."
!else
	MessageBox MB_OK "An existing version of the libjpeg-turbo SDK for Visual C++ 64-bit or the TurboJPEG SDK is already installed.  Please uninstall it first."
!endif
	quit

	notexists:
	SetOutPath $SYSDIR
!ifdef GCC
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\libturbojpeg.dll"
!else
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\${BUILDDIR}turbojpeg.dll"
!endif
	SetOutPath $INSTDIR\bin
!ifdef GCC
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\libturbojpeg.dll"
!else
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\${BUILDDIR}turbojpeg.dll"
!endif
!ifdef GCC
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\libjpeg-62.dll"
!else
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\${BUILDDIR}jpeg62.dll"
!endif
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\${BUILDDIR}cjpeg.exe"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\${BUILDDIR}djpeg.exe"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\${BUILDDIR}jpegtran.exe"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\${BUILDDIR}tjbench.exe"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\${BUILDDIR}rdjpgcom.exe"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\${BUILDDIR}wrjpgcom.exe"
	SetOutPath $INSTDIR\lib
!ifdef GCC
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\libturbojpeg.dll.a"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\libturbojpeg.a"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\libjpeg.dll.a"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\libjpeg.a"
!else
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\${BUILDDIR}turbojpeg.lib"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\${BUILDDIR}turbojpeg-static.lib"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\${BUILDDIR}jpeg.lib"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\${BUILDDIR}jpeg-static.lib"
!endif
	SetOutPath $INSTDIR\lib\pkgconfig
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\pkgscripts\libjpeg.pc"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\pkgscripts\libturbojpeg.pc"
	SetOutPath $INSTDIR\lib\cmake\libjpeg-turbo
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\pkgscripts\libjpeg-turboConfig.cmake"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\pkgscripts\libjpeg-turboConfigVersion.cmake"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\win\libjpeg-turboTargets.cmake"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\win\libjpeg-turboTargets-release.cmake"
!ifdef JAVA
	SetOutPath $INSTDIR\classes
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\java\turbojpeg.jar"
!endif
	SetOutPath $INSTDIR\include
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\jconfig.h"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\jerror.h"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\jmorecfg.h"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\jpeglib.h"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\turbojpeg.h"
	SetOutPath $INSTDIR\doc
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\README.ijg"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\README.md"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\LICENSE.md"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\example.txt"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\libjpeg.txt"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\structure.txt"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\usage.txt"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\wizard.txt"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\tjexample.c"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\java\TJExample.java"
!ifdef GCC
	SetOutPath $INSTDIR\man\man1
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\cjpeg.1"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\djpeg.1"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\jpegtran.1"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\rdjpgcom.1"
	File "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo\wrjpgcom.1"
!endif

	WriteRegStr HKLM "SOFTWARE\libjpeg-turbo64 2.1.1" "Install_Dir" "$INSTDIR"

	WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\libjpeg-turbo64 2.1.1" "DisplayName" "libjpeg-turbo SDK v2.1.1 for Visual C++ 64-bit"
	WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\libjpeg-turbo64 2.1.1" "UninstallString" '"$INSTDIR\uninstall_2.1.1.exe"'
	WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\libjpeg-turbo64 2.1.1" "NoModify" 1
	WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\libjpeg-turbo64 2.1.1" "NoRepair" 1
	WriteUninstaller "uninstall_2.1.1.exe"
SectionEnd

Section "Uninstall"
!ifdef WIN64
	${If} ${RunningX64}
	${DisableX64FSRedirection}
	${Endif}
!endif

	SetShellVarContext all

	DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\libjpeg-turbo64 2.1.1"
	DeleteRegKey HKLM "SOFTWARE\libjpeg-turbo64 2.1.1"

!ifdef GCC
	Delete $INSTDIR\bin\libjpeg-62.dll
	Delete $INSTDIR\bin\libturbojpeg.dll
	Delete $SYSDIR\libturbojpeg.dll
	Delete $INSTDIR\lib\libturbojpeg.dll.a
	Delete $INSTDIR\lib\libturbojpeg.a
	Delete $INSTDIR\lib\libjpeg.dll.a
	Delete $INSTDIR\lib\libjpeg.a
!else
	Delete $INSTDIR\bin\jpeg62.dll
	Delete $INSTDIR\bin\turbojpeg.dll
	Delete $SYSDIR\turbojpeg.dll
	Delete $INSTDIR\lib\jpeg.lib
	Delete $INSTDIR\lib\jpeg-static.lib
	Delete $INSTDIR\lib\turbojpeg.lib
	Delete $INSTDIR\lib\turbojpeg-static.lib
!endif
	Delete $INSTDIR\lib\pkgconfig\libjpeg.pc
	Delete $INSTDIR\lib\pkgconfig\libturbojpeg.pc
	Delete $INSTDIR\lib\cmake\libjpeg-turbo\libjpeg-turboConfig.cmake
	Delete $INSTDIR\lib\cmake\libjpeg-turbo\libjpeg-turboConfigVersion.cmake
	Delete $INSTDIR\lib\cmake\libjpeg-turbo\libjpeg-turboTargets.cmake
	Delete $INSTDIR\lib\cmake\libjpeg-turbo\libjpeg-turboTargets-release.cmake
!ifdef JAVA
	Delete $INSTDIR\classes\turbojpeg.jar
!endif
	Delete $INSTDIR\bin\cjpeg.exe
	Delete $INSTDIR\bin\djpeg.exe
	Delete $INSTDIR\bin\jpegtran.exe
	Delete $INSTDIR\bin\tjbench.exe
	Delete $INSTDIR\bin\rdjpgcom.exe
	Delete $INSTDIR\bin\wrjpgcom.exe
	Delete $INSTDIR\include\jconfig.h
	Delete $INSTDIR\include\jerror.h
	Delete $INSTDIR\include\jmorecfg.h
	Delete $INSTDIR\include\jpeglib.h
	Delete $INSTDIR\include\turbojpeg.h
	Delete $INSTDIR\uninstall_2.1.1.exe
	Delete $INSTDIR\doc\README.ijg
	Delete $INSTDIR\doc\README.md
	Delete $INSTDIR\doc\LICENSE.md
	Delete $INSTDIR\doc\example.txt
	Delete $INSTDIR\doc\libjpeg.txt
	Delete $INSTDIR\doc\structure.txt
	Delete $INSTDIR\doc\usage.txt
	Delete $INSTDIR\doc\wizard.txt
	Delete $INSTDIR\doc\tjexample.c
	Delete $INSTDIR\doc\TJExample.java
!ifdef GCC
	Delete $INSTDIR\man\man1\cjpeg.1
	Delete $INSTDIR\man\man1\djpeg.1
	Delete $INSTDIR\man\man1\jpegtran.1
	Delete $INSTDIR\man\man1\rdjpgcom.1
	Delete $INSTDIR\man\man1\wrjpgcom.1
!endif

	RMDir "$INSTDIR\include"
	RMDir "$INSTDIR\lib\pkgconfig"
	RMDir "$INSTDIR\lib\cmake\libjpeg-turbo"
	RMDir "$INSTDIR\lib\cmake"
	RMDir "$INSTDIR\lib"
	RMDir "$INSTDIR\doc"
!ifdef GCC
	RMDir "$INSTDIR\man\man1"
	RMDir "$INSTDIR\man"
!endif
!ifdef JAVA
	RMDir "$INSTDIR\classes"
!endif
	RMDir "$INSTDIR\bin"
	RMDir "$INSTDIR"

SectionEnd
