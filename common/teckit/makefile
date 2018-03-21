#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------

#--------------------------------------
# TARGET
#--------------------------------------
TARGET = main

#--------------------------------------
# Source files
#--------------------------------------
FILE_MAP  = punctuation.tec   \
            ascii-uipa.tec    \
            test.tec          \
            asciitogreek.tec  \
            Preeti.tec        \
            zht2zhs.tec

#--------------------------------------
# Paths
#--------------------------------------
 PATH_TECKIT = c:\p\teckit\32-bit

#--------------------------------------
# Tools
#--------------------------------------
 PRG_COMPILE = $(PATH_TECKIT)\teckit_compile.exe

#--------------------------------------
# Build Control
#--------------------------------------
$(TARGET).log: $(FILE_MAP)
   dir /od >> $(TARGET).log

ascii-uipa.tec: ascii-uipa.map $(PRG_COMPILE)
  $(PRG_COMPILE) ascii-uipa.map

asciitogreek.tec: asciitogreek.map $(PRG_COMPILE)
  $(PRG_COMPILE) asciitogreek.map

cjk-pinyin.tec: cjk-pinyin.map $(PRG_COMPILE)
  $(PRG_COMPILE) cjk-pinyin.map

punctuation.tec: punctuation.map $(PRG_COMPILE)
  $(PRG_COMPILE) punctuation.map

Preeti.tec: Preeti.map $(PRG_COMPILE)
  $(PRG_COMPILE) Preeti.map

zht2zhs.tec: zht2zhs.map $(PRG_COMPILE)
  $(PRG_COMPILE) zht2zhs.map

test.tec: test.map $(PRG_COMPILE)
  $(PRG_COMPILE) test.map

