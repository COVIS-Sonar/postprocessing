help:
	@echo "make download_all				Downloads all of the COVIS test data.  Requires wget"

TEST_FILES = COVIS-20191024T003346-diffuse3.7z

download_all: $(TEST_FILES)

COVIS-*:
	wget $@

.phony:  help download_all
