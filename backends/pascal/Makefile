# Makefile -- kh-sim Pascal backend (KH-007), phase-7f canonical-physics targets.
# Companion to build.ps1 (the canonical Windows build tool). This
# Makefile provides the `make test-canonical` idiom used by CHILD-pascal.md S1
# for POSIX/native FPC parity.
#
# Requirements: Free Pascal Compiler (fpc) in PATH.
#
# Targets:
#   make server          build the HTTP server (src/kh_server.pas -> build/)
#   make canonical       build the phase-7f test (tests/test_canonical_physics.pas)
#   make test-canonical  build + run the phase-7f test (server must be running)
#   make validate        build the phase-7d validation binary
#   make clean

FPC      ?= fpc
# Reference build flags (no perf flags): range + overflow checks, large stack.
# Mirrors CHILD-pascal.md S1: -O- -Cr -Co -Cs16777216.
FPCFLAGS ?= -O- -Cr -Co -Cs16777216 -Fusrc -Fubuild -FUbuild -FEbuild
BUILD    := build

$(BUILD):
	mkdir -p $(BUILD)

server: $(BUILD)
	$(FPC) $(FPCFLAGS) src/kh_server.pas

validate: $(BUILD)
	$(FPC) $(FPCFLAGS) tests/kh_validate.pas

canonical: $(BUILD)
	$(FPC) $(FPCFLAGS) -Futests tests/test_canonical_physics.pas

# Run the canonical-physics testbed. The HTTP server must already be listening
# on port 8005 (start it with: ./build/kh_server &).
test-canonical: canonical
	./$(BUILD)/test_canonical_physics

clean:
	rm -f $(BUILD)/*.o $(BUILD)/*.ppu $(BUILD)/test_canonical_physics \
	      $(BUILD)/kh_server $(BUILD)/kh_validate

.PHONY: server validate canonical test-canonical clean
