help:
	@echo "make or make download 			Downloads all of the COVIS test data.  Requires wget"

TEST_FILES = COVIS-20191024T003346-diffuse3.7z

default: download

download: $(TEST_FILES)

COVIS-*:
	wget $@

.phony:  help download
