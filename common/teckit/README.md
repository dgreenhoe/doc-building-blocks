#============================================================================
# Daniel J. Greenhoe
# Readme file for 
# https://github.com/dgreenhoe/teckitmaps/
#============================================================================
#---------------------------------------
# What is a TECkit map?
#---------------------------------------
A TECkit map is a Unicode mapping file which, after compiling, can be used with 
the XeLaTeX software to map one or more Unicode points to another Unicode point.
  scripts.sil.org/teckit
Note that before being used, a TECkit .map file must be compiled to a .tec file.

#---------------------------------------
# TECkit maps in this repository
#---------------------------------------
  * punctuation.map: punctuation and PinYin (for Chinese phonetics) mappings
  * ascii-uipa.map:  IPA symbols
  * zht2zhs.map:     Traditional Chinese character to Chinese Simplified character mapping

#---------------------------------------
# Example 1 [the emdash]
#---------------------------------------
In an .tex file, one might well want an "em dash".
An "em dash" is simply a long dash that is approximately 1em wide.
"1 em wide" means that it is as wide as a capital Latin letter (ABC...)
might be expected to be.
In the Unicode standard, the em dash is code point U+2014:
  www.fileformat.info/info/Unicode/char/2014/
  http://www.Unicode.org/charts/PDF/U2000.pdf
In a .tex source file, entering the em dash directly is arguably inconvenient.
One solution to this problem is to define an entry in a TECkit map file as follows:
  hyphen_minus hyphen_minus hyphen_minus <> em_dash  ;---<> U+2014
Note that the ";" and everything to the right is a comment.
Now, whenever using a font that specifies the TECkit map, 
a "---" in the .tex source will be typeset as an em dash (U+2014).
In this repository, such mapping is performed by the TECkit map "punctuation".
To specify use of this TECkit map, one can do this in the Unicode font setup:
 \defaultfontfeatures{%
   Mapping = punctuation,
   }

#---------------------------------------
# Example 2 [PinYin for Chinese phonetics]
#---------------------------------------
Spoken Mandarin Chinese has about 401 possible syllable sounds ("phones").
This isn't many. But each syllable can be one of five tones.
These tones are referred to as "first tone", "second tone", "third tone",
"fourth tone", and "neutral":
  https://www.researchgate.net/project/Language-Toolbox
The language can be described using any of a number of phonetic systems.
One of the most popular is PinYin.
In the PinYin system, 
  * first tone  is marked with a horizontal dash over one of the letters in the syllable, 
  * second tone is marked with a rising dash,
  * third tone  is marked with a "v" shaped dash
  * fourth tone is marked with a falling dash.
For convenient typesetting of such tones into Unicode points, one could do this:
  FSL HYP            latin_small_letter_a <> latin_small_letter_a_with_macron     ;/-a <> U+0101
  FSL apostrophe     latin_small_letter_a <> latin_small_letter_a_with_acute      ;/'a <> U+00E1
  FSL less_than_sign latin_small_letter_a <> latin_small_letter_a_with_circumflex ;/<a <> U+00E2
  FSL grave_accent   latin_small_letter_a <> latin_small_letter_a_with_grave      ;/`a <> U+00E0
Then, in the .tex source, one can type "sh/-an", which will be typeset as 
"shan" with a bar over the top of the "a", which is a way to say "mountain" in 
Mandarin Chinese.
In this repository, such mapping is performed by the TECkit map
  punctuation.map

#---------------------------------------
# Example 3 [IPA phonetics]
#---------------------------------------
The International Phonetic Alphabet (IPA) contains many special symbols that can be used for 
the phonetic description of many many languages. 
  http://unicode.org/charts/PDF/U0250.pdf
These special symbols are contained in some fonts such as 
  "GnuFree Sans-Serif"  (http://www.gnu.org/s/freefont/)
One way to conveniently typeset these glyphs is by using a TECkit mapping file with entries such as
  Define FSL solidus ; U+002F ("/", "forward slash")
  FSL latin_small_letter_t latin_capital_letter_S   <> latin_small_letter_tesh_digraph  ; /tS  <> U+02A7  
Then, entering in a .tex file the sequence "/tS" will be rendered with 
the single Unicode glyph U+02A7,
which looks like a "t" and an "S" squished together into a single glyph and 
represents the "voiceless postalveolar affricate" sound (e.g. the "ch" sound in English).
In this repository, such mappings are provided in "ascii-uipa.map".
An appropriate font can be assigned the mappings using something like this:
  \newcommand{\fntipa} {\fntFreeSans\addfontfeatures{Mapping=/dan/r/common/teckit/ascii-uipa}}

#---------------------------------------
# Example 4 [Traditional to Simplified Chinese characters]
#---------------------------------------
Currently both traditional and simplified Chinese characters are 
both used by large numbers of people in the world.
Therefore, it is "best" that documents typeset "in Chinese" 
have both a traditional and simplified versions available.
One solution to this problem is to typeset the document using traditional characters,
and then use a TECkit map to map the appropriate characters into the simplified characters.
In this repository, such a mapping is provided by "zht2zhs.map".
