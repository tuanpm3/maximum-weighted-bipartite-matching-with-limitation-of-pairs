.PHONY: clean All

All: release
debug:
	@echo "----------Building project:[ hungarian - Debug ]----------"
	@"$(MAKE)" -f  "hungarian_debug.mk"
release:
	@echo "----------Building project:[ hungarian - Release ]----------"
	@"$(MAKE)" -f  "hungarian_release.mk"
clean:
	@echo "----------Cleaning project:[ hungarian - Release ]----------"
	@"$(MAKE)" -f  "hungarian_release.mk" clean
	@"$(MAKE)" -f  "hungarian_debug.mk" clean
