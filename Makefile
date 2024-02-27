include Makefile.config

main: 
	$(MAKE) -C src

testAll:
	$(MAKE) run -C tests

clean:
	@echo "Removing all object files"
	@ -find . -name "*.o" -exec rm {} \;
	@echo "Removing $(BUILD_PATH) dirs"
	@ -rm -rf src/$(BUILD_PATH)
	@ -rm -rf tests/$(BUILD_PATH)
