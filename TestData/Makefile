help:
	@echo "make download 			Downloads all of the COVIS test data.  Requires wget"

TEST_FILES = COVIS-20191024T003346-diffuse3.tar.gz

download: $(TEST_FILES)

%.tar.gz:
	wget https://covis-raw.s3.wasabisys.com/2019/10/24/$@

.phony:  help download
