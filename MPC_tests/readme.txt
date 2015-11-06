MPC tests

Generate JSim mod files and diff them with mod files in ref subdirectory 
to confirm no regressions in new code:

1. Generate JSim mod files:
java -jar MPC.jar hux_module.mpc
java -jar MPC.jar huxley_name_contains.mpc
java -jar MPC.jar Ado2Ino3comp.mpc
java -jar MPC.jar example.mpc

2. Diff output with reference files:
 diff hux_module.mod ref/hux_module.mod
 diff huxley_name_contains.mod ref/huxley_name_contains.mod
 diff Ado2Ino3comp.mod ref/Ado2Ino3comp.mod
 diff example.mod ref/example.mod

