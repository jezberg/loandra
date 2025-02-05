default: release-loginfo

define update_compile_commands =
@rm -f compile_commands.json
@ln -s $(PWD)/build/$(subst libonly-,,$@)/compile_commands.json compile_commands.json
endef

release-loginfo \
release-loginfo-noffa \
release-loginfo-cadical \
release-logverbose \
release-logtrace \
debug-loginfo \
debug-logverbose \
debug-logtrace \
release-profiling \
static-release-logverbose \
:
	cmake -S . --preset=$@
	$(update_compile_commands)
	cmake --build build/$@

libonly-debug-logverbose \
libonly-release-nologging \
:
	cmake -S . --preset=$(subst libonly-,,$@)
	$(update_compile_commands)
	cmake --build  --preset=$@

clean::
	@rm -rf build compile_commands.json

.PHONY: \
default \
release-loginfo \
release-loginfo-noffa \
release-loginfo-cadical \
release-logverbose \
release-logtrace \
debug-loginfo \
debug-logverbose \
debug-logtrace \
release-profiling \
static-release-logverbose \
libonly-debug-logverbose \
libonly-release-nologging \
clean
