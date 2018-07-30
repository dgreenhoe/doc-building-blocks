#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------

#--------------------------------------
# Files
#--------------------------------------


#--------------------------------------
# Build Control
#--------------------------------------
#--------------------------------------
# Commands
#--------------------------------------
all:
  cd bsplines
    make
  cd ..\dsp
    make
  cd ..\xwsd
    make
#  cd ..\xpsd
#    make
  cd ..\xsams
    make
  cd ..\estimation
    make
  cd ..\dcsd
    make
  cd ..\A290826
    make
  cd ..\lmc
    make
  cd ..\elife
    make
  cd ..\zsounds
    make
  cd ..\xsro
    make
  cd ..\trigsys
    make
  cd ..\hasse
    make
  cd ..\A289358
    make
  cd ..\2015ssp
    make
  cd ..\abcssp
    make
  cd ..\2016apcssp
    make
  cd ..\nil
    make
  cd ..\2015larc
    make
  cd ..\2015pds
    make
  cd ..\2014flat
    make
  cd ..\2014sto
    make
  cd ..\2014punity
    make
  cd ..\2014td
    make
  cd ..\2014plat
   make
  cd ..

cleanall:
  cd frames
    make clean
  cd ..\xwsd
    make clean
  cd ..\xpsd
    make clean
  cd ..\xsams
    make clean
  cd ..\dcsd
    make clean
  cd ..\bsplines
    make clean
  cd ..\estimation
    make clean
  cd ..\A290826
    make clean
  cd ..\xsro
    make clean
  cd ..\lmc
    make clean
  cd ..\elife
    make clean
  cd ..\zsounds
    make clean
  cd ..\trigsys
    make clean
  cd ..\hasse
    make clean
  cd ..\A289358
    make clean
  cd ..\2015ssp
    make clean
  cd ..\abcssp
    make clean
  cd ..\2016apcssp
    make clean
  cd ..\nil
    make clean
  cd ..\2015larc
    make clean
  cd ..\2015pds
    make clean
  cd ..\2014flat
    make clean
  cd ..\2014sto
    make clean
  cd ..\2014punity
    make clean
  cd ..\2014td
    make clean
  cd ..\2014plat
    make clean
  cd ..

new:
   make clean
   make all

print:
  c:\p\gs902\bin\gswin32c -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$(TARGET)xsd.pdf ..\2014flat\2014flat.pdf ..\2014plat\2014plat.pdf ..\2014td\2014td.pdf ..\2014punity\2014punity.pdf 

ap:
  cd apengzhtA5
    make
  cd ..\apengzhtA4
    make
  cd ..\apengzhsA5
    make
  cd ..\apengzhsA4
    make
  cd ..\apengA5
    make
  cd ..\apengA4
    make
  cd ..\apzhsA5
    make
  cd ..\apzhsA4
    make
  cd ..\apzhtA5
    make
  cd ..\apzhtA4
    make
  cd ..\apzht46
    make
  cd ..\apzhs46
    make
  cd ..\apeng46
    make
  cd ..


backup:
  zip -o -9 -r xsd.zip * -x *\plots\*.tex *\bingo*.tex *\old\* *\wavr_old\* tmp*.* t.* *.mp4 *.ppm *.xcf *.odp *.odt *.ttf *.mp3 *.wav *.synctex *.ctx *.zip *.pdf *.ps *.dvi *.aux *.eps *.exe *.prn *.obj *.tmp *.bak *.bat *.tds *.png *.jpg *.jpeg *.gif *.bmp *.log *.js *.css *.swf *.xdv *.dat *.tec *.otf *.*dx *.url *.out *.tar *.gz *.0 *.txt *.htm *.html *\*files\* *\mainx.* *\s\* *\trash\* trash\* *\save\* *\save2\* *\deprecated\* deprecated\* tmp\* ec\*.tex *.xlg *.ind *.ilg

archiveap:
  #copy /v apengzhtA5\apengzhtA5.pdf  apengzhtA5\pdfs\AncientPath_engzhtA5_20180116x.pdf
  #copy /v apengzhtA4\apengzhtA4.pdf  apengzhtA5\pdfs\AncientPath_engzhtA4_20180116x.pdf
  #copy /v apengzhsA5\apengzhsA5.pdf  apengzhtA5\pdfs\AncientPath_engzhsA5_20180116x.pdf
  #copy /v apengzhsA4\apengzhsA4.pdf  apengzhtA5\pdfs\AncientPath_engzhsA4_20180116x.pdf
  #copy /v apengA5\apengA5.pdf        apengzhtA5\pdfs\AncientPath_engA5_20180116x.pdf   
  #copy /v apengA4\apengA4.pdf        apengzhtA5\pdfs\AncientPath_engA4_20180116x.pdf   
  #copy /v apzhsA5\apzhsA5.pdf        apengzhtA5\pdfs\AncientPath_zhsA5_20180116x.pdf   
  #copy /v apzhsA4\apzhsA4.pdf        apengzhtA5\pdfs\AncientPath_zhsA4_20180116x.pdf   
  #copy /v apzhtA5\apzhtA5.pdf        apengzhtA5\pdfs\AncientPath_zhtA5_20180116x.pdf   
  #copy /v apzhtA4\apzhtA4.pdf        apengzhtA5\pdfs\AncientPath_zhtA4_20180116x.pdf   
  #copy /v apzht46\apzht46.pdf         apengzhtA5\pdfs\AncientPath_zht46_20180116x.pdf   
  #copy /v apzhs46\apzhs46.pdf         apengzhtA5\pdfs\AncientPath_zhs46_20180116x.pdf   
  #copy /v apeng46\apeng46.pdf         apengzhtA5\pdfs\AncientPath_eng46_20180116x.pdf   
  #copy /v apengA5\apengA5.pdf         apengzhtA5\pdfs\AncientPath_engA5_20180116x.pdf   
  #copy /v apengA4\apengA4.pdf         apengzhtA5\pdfs\AncientPath_engA4_20180116x.pdf   
  #
  #copy /v apeng46\apeng46.pdf        apengzhtA5\pdfs\AncientPath_eng46-010_20180116x.pdf   
  #copy /v apengA5\apengA5.pdf        apengzhtA5\pdfs\AncientPath_engA5-010_20180116x.pdf   
  #copy /v apzhtA5\apzhtA5.pdf        apengzhtA5\pdfs\AncientPath_zhtA5-010_20180116x.pdf   
  #copy /v apzhsA5\apzhsA5.pdf        apengzhtA5\pdfs\AncientPath_zhsA5-011_20180117.pdf   
  #copy /v apengzhtA5\apengzhtA5.pdf  apengzhtA5\pdfs\AncientPath_engzhtA5-011_20180117.pdf
  #copy /v apengzhsA5\apengzhsA5.pdf  apengzhtA5\pdfs\AncientPath_engzhsA5-011_20180117.pdf
  #copy /v apeng46w\apeng46.pdf       apengzhtA5\pdfs\AncientPath_eng46w-010_20180116x.pdf   
  #
  #copy /v apzhsA5\apzhsA5.pdf       ap\pdfs\AncientPath_zhsA5-012_20180117.pdf
  #copy /v apzhtA5\apzhtA5.pdf       ap\pdfs\AncientPath_zhtA5-012_20180117.pdf
  #copy /v apengzhtA5\apengzhtA5.pdf ap\pdfs\AncientPath_engzhtA5-012_20180117.pdf
  #copy /v apengzhsA5\apengzhsA5.pdf ap\pdfs\AncientPath_engzhsA5-012_20180117.pdf
  #copy /v apengzhtA4\apengzhtA4.pdf ap\pdfs\AncientPath_engzhtA4-012_20180117.pdf
  #copy /v apengzhsA4\apengzhsA4.pdf ap\pdfs\AncientPath_engzhsA4-012_20180117.pdf
  #copy /v apzhsA4\apzhsA4.pdf       ap\pdfs\AncientPath_zhsA4-012_20180117.pdf
  #copy /v apzhtA4\apzhtA4.pdf       ap\pdfs\AncientPath_zhtA4-012_20180117.pdf
  #copy /v apzhs46\apzhs46.pdf       ap\pdfs\AncientPath_zhs46-012_20180117.pdf
  #copy /v apzht46\apzht46.pdf       ap\pdfs\AncientPath_zht46-012_20180117.pdf
  #copy /v apeng46\apeng46.pdf       ap\pdfs\AncientPath_eng46-012_20180118.pdf
  #
  copy /v apzhsA5\apzhsA5.pdf       ap\pdfs\AncientPath_zhsA5-014_20180118.pdf
  copy /v apzhtA5\apzhtA5.pdf       ap\pdfs\AncientPath_zhtA5-014_20180118.pdf
  copy /v apengzhtA5\apengzhtA5.pdf ap\pdfs\AncientPath_engzhtA5-014_20180118.pdf
  copy /v apengzhsA5\apengzhsA5.pdf ap\pdfs\AncientPath_engzhsA5-014_20180118.pdf
  copy /v c:\dan\r\common\graphics\watercraft\4char.pdf ap\pdfs\
  copy /v c:\dan\r\common\graphics\watercraft\4char.jpg ap\pdfs\



