#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------
#--------------------------------------
# TARGET
#--------------------------------------
TARGET = tl

#--------------------------------------
# Files
#--------------------------------------
FILES_OBJ = main.o command.o estimate_amplitude.o estimate_induction.o help.o sequencer.o traffic.o
FILES_H   = main.h command.h estimate_amplitude.h estimate_induction.h help.h sequencer.h traffic.h

#--------------------------------------
# directories
#--------------------------------------
#DIR_BIN = c:\p\bcc\bin      # binaries (execuatables) directory
#DIR_INC = c:\p\bcc\Include  # include files directory
#DIR_LIB = c:\p\bcc\Lib      # libraries directory

#--------------------------------------
# Programs
#--------------------------------------
PRG_COMPILE = g++
PRG_LINK    = g++

#--------------------------------------
# master build control
#--------------------------------------
$(TARGET).exe: $(FILES_OBJ) $(FILES_H) makefile
  time /T
  $(PRG_LINK)  -o $(TARGET).exe $(FILES_OBJ)
  dir $(TARGET).exe
#  @echo Done ... sto/graphics pdfs are ready!

#--------------------------------------
# implicit build control
#--------------------------------------
.cpp.o:
	$(PRG_COMPILE) -c $&.cpp

#--------------------------------------
# commands
#--------------------------------------
clean:
  del *.o

scrub:
  make clean
  del *.exe

new:
  del *.obj
  make

zip:
  zip -o -9 -r $(TARGET).zip * -x *.obj *.tds *.exe *.dat *.sty *.zip *.bat *.pdf tmp\*
  unzip -l $(TARGET).zip
  dir *.zip
