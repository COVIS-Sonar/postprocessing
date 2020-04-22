help:
	@echo "make download 			Downloads all of the COVIS test data.  Requires wget"

TEST_FILES = 2019/10/24/COVIS-20191024T003346-diffuse3.tar.gz \
							2019/10/24/COVIS-20191024T000002-imaging1.tar.gz \
							2020/04/21/COVIS-20200421T133002-diffuse1.tar.gz \
							2019/10/24/COVIS-20191024T003346-diffuse3.7z \
							2019/10/24/COVIS-20191024T000002-imaging1.7z \
							2020/04/21/COVIS-20200421T133002-diffuse1.7z


download: $(TEST_FILES)

%.tar.gz:
	mkdir -p $(dir $@)
	cd $(dir $@) && wget https://covis-raw.s3.wasabisys.com/$@

## Since the 7z files aren't available on the public internet,
## make them locally instead.
## Computationally expensive
%.7z: %.tar.gz
	mkdir -p $(dir $@)
	$(eval TMP := $(shell mktemp -d))
	cd $(TMP) && tar xzvf $(abspath $^) | xargs 7z a $(abspath $@)
	rm -rf $(TMP)


.phony:  help download
