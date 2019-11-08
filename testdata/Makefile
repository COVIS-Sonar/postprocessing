help:
	@echo "make download 			Downloads all of the COVIS test data.  Requires wget"

TEST_FILES = 2019/10/24/COVIS-20191024T003346-diffuse3.tar.gz \
							2019/10/24/COVIS-20191024T000002-imaging1.tar.gz \
							2019/10/24/COVIS-20191024T003346-diffuse3.7z \
							2019/10/24/COVIS-20191024T000002-imaging1.7z


download: $(TEST_FILES)

%.tar.gz:
	mkdir -p $(dir $@)
	cd $(dir $@) && wget https://covis-raw.s3.wasabisys.com/$@

%.7z: %.tar.gz
	mkdir -p $(dir $@)
	tar xz $^ -O | 7z a -si $@


.phony:  help download
