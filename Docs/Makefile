
COVIS_DATA_FORMAT=CovisDataFormat.md

all: $(COVIS_DATA_FORMAT:.md=.txt)

$(COVIS_DATA_FORMAT:.md=.txt): $(COVIS_DATA_FORMAT)
	pandoc $^ -t plain -o $@


.PHONY: all
