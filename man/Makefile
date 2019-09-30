#!/bin/make -f

TXT  := RELEASE_NOTES
## 1   Executable programs or shell commands
TXT1 := dftd4 cosmosolv cosmothermrd thermo xtb
## 2   System calls
TXT2 :=
## 3   Library calls (functions within program libraries)
TXT3 := thermo distort
## 4   Special files (usually found in /dev)
TXT4 :=
## 5   File formats and conventions eg /etc/passwd
TXT5 :=
## 6   Games
TXT6 :=
## 7   Miscellaneous (including macro packages and conventions)
## here we place all benchmarks in our arsenal
TXT7 := xcontrol
#MAN7 += ACONF ADIM6 AHB21 AL2X6 ALK8 ALKBDE10 Amino20x4
#MAN7 += BH76 BH76RC BHDIV10 BHPERI BHROT27 BSR36 BUT14DIOL
#MAN7 += Bühl-3d Bühl-4d Bühl-5d
#MAN7 += C15 C60ISO CARBHB12 CDIE20 CHB6
#MAN7 += DARC DC13 DIPCS10
#MAN7 += FH51
#MAN7 += G21EA G21IP G2RC GrubbsII
#MAN7 += HAL59 HEAVY28 HEAVYSB11 HMGB11
#MAN7 += ICONF IDISP IL16 IL-2013 INV24 ISO34 ISOL24
#MAN7 += LMGB35
#MAN7 += MB16-43 MCONF MOR41
#MAN7 += NBPRC NCIBLIND10
#MAN7 += PA26 PArel PCONF21 PNICO23 PX13
#MAN7 += RC21 RG18 RSE43 R160x6
#MAN7 += S22 S22x5 S66 S66x8 S30L SCONF SIE4x4
#MAN7 += TAUT15 TOMC6
#MAN7 += UPU23
#MAN7 += W4-11 WATER27 WCPT18 WCCR10
#MAN7 += YBDE18
#MAN7 += X40x10

MAN1 := $(patsubst %,man1/%.1,$(TXT1))
MAN2 := $(patsubst %,man2/%.2,$(TXT2))
MAN3 := $(patsubst %,man3/%.3,$(TXT3))
MAN4 := $(patsubst %,man4/%.4,$(TXT4))
MAN5 := $(patsubst %,man5/%.5,$(TXT5))
MAN6 := $(patsubst %,man6/%.6,$(TXT6))
MAN7 := $(patsubst %,man7/%.7,$(TXT7))
MAN := $(MAN1)
MAN += $(MAN2)
MAN += $(MAN3)
MAN += $(MAN4)
MAN += $(MAN5)
MAN += $(MAN6)
MAN += $(MAN7)

PDF  := $(patsubst %,%.1.pdf,$(TXT1))
PDF  += $(patsubst %,%.2.pdf,$(TXT2))
PDF  += $(patsubst %,%.3.pdf,$(TXT3))
PDF  += $(patsubst %,%.4.pdf,$(TXT4))
PDF  += $(patsubst %,%.5.pdf,$(TXT5))
PDF  += $(patsubst %,%.6.pdf,$(TXT6))
PDF  += $(patsubst %,%.7.pdf,$(TXT7))

HTML := $(patsubst %,%.1.html,$(TXT1))
HTML += $(patsubst %,%.2.html,$(TXT2))
HTML += $(patsubst %,%.3.html,$(TXT3))
HTML += $(patsubst %,%.4.html,$(TXT4))
HTML += $(patsubst %,%.5.html,$(TXT5))
HTML += $(patsubst %,%.6.html,$(TXT6))
HTML += $(patsubst %,%.7.html,$(TXT7))

vpath %.txt .

.PHONY: all

all: $(MAN) $(PDF) $(HTML)

ASCIIDOC := asciidoc
A2X := a2x
ASCIIFLAGS := --doctype manpage --format manpage

DBLATEX  := dblatex
DBLATEXFLAGS := --style simple

$(HTML): %.html: %.txt
	$(ASCIIDOC) $(ASCIIDOCFLAGS) --backend html $<

$(PDF): %.pdf: %.xml
	$(DBLATEX) $(DBLATEXFLAGS) $<

%.xml: %.txt
	$(ASCIIDOC) $(ASCIIDOCFLAGS) --backend docbook $<

$(MAN1): man1/%.1: %.1.txt
	$(A2X) $(ASCIIFLAGS) -D man1 $<

$(MAN2): man2/%.2: %.2.txt
	$(A2X) $(ASCIIFLAGS) -D man2 $<

$(MAN3): man3/%.3: %.3.txt
	$(A2X) $(ASCIIFLAGS) -D man3 $<

$(MAN4): man4/%.4: %.4.txt
	$(A2X) $(ASCIIFLAGS) -D man4 $<

$(MAN5): man5/%.5: %.5.txt
	$(A2X) $(ASCIIFLAGS) -D man5 $<

$(MAN6): man6/%.6: %.6.txt
	$(A2X) $(ASCIIFLAGS) -D man6 $<

$(MAN7): man7/%.7: %.7.txt
	$(A2X) $(ASCIIFLAGS) -D man7 $<

.PHONY: clean

clean:
	$(RM) $(MAN1) $(MAN7)
