::===========================================================================
:: Daniel J. Greenhoe
:: batch file for converting images to 10k-15k size
:: uses ImageMagik "convert.exe"
:: see http://www.imagemagick.org/
::===========================================================================

convert -density 72 -resize 600x600 -quality 90 FirstThanksgiving_wiki.jpg  small\FirstThanksgiving_wiki.jpg
dir tiny /od 
dir small /od
