all: build/main.pdf

# hier Python-Skripte:
build/emission.pdf: emission_neu.py matplotlibrc header-matplotlib.tex | build
	TEXINPUTS=$$(pwd): python emission_neu.py

build/Zink.pdf: absorption_neu.py matplotlibrc header-matplotlib.tex | build
	TEXINPUTS=$$(pwd): python absorption_neu.py

build/bra.pdf: bragg.py matplotlibrc header-matplotlib.tex | build
	TEXINPUTS=$$(pwd): python bragg.py

build/plot.pdf: vorb_v602.py matplotlibrc header-matplotlib.tex | build
	TEXINPUTS=$$(pwd): python vorb_v602.py

# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/main.pdf: build/emission.pdf build/bra.pdf build/Zink.pdf build/plot.pdf

build/main.pdf: FORCE | build
	  TEXINPUTS=build: \
	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	main.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
