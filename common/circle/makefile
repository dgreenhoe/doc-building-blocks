#-----------------------------------------------------------------------------
# Project Makefile
# Daniel Greenhoe
#-----------------------------------------------------------------------------

#--------------------------------------
# Directories
#--------------------------------------
DIR_COMP  = e:\bcc55\bin
DIR_LINK  = e:\bcc55\bin
DIR_INCL  = e:\bcc55\include
DIR_LIB   = e:\bcc55\lib

#--------------------------------------
# Tools
#--------------------------------------
COMPILER  = $(DIR_COMP)\bcc32
LINKER    = $(DIR_LINK)\bcc32
CFLAGS    = -c -K -O2 -I$(DIR_INCL)
LFLAGS    = -M -L$(DIR_LIB)

#--------------------------------------
# Files
#--------------------------------------
TARGET    = drawc
FILE_EXE  = $(TARGET).exe
FILE_MAP  = $(TARGET).map
FILE_OBJ  =              \
            drawc.obj     

#--------------------------------------
# Build Control
#--------------------------------------
$(FILE_EXE): $(FILE_OBJ)
   $(LINKER) -e$(TARGET) $(LFLAGS) $(FILE_OBJ)

.c.obj:
   $(COMPILER) $(CFLAGS) {$< }

#--------------------------------------
# Commands
#--------------------------------------
clean:
   del *.obj
   del *.tds
   del *.map

new:
   make clean
   make

backup:
   make clean
   pkzip -z -x*.dat -xc1 -xc2 $(TARGET).zip *.*
   
go:
   make
   $(TARGET) -rdx 200 36 -q

