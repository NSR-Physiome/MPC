MPC tests

Generate JSim mod files and diff them with mod files in ref subdirectory 
to confirm no regressions in new code:

1. Generate JSim mod files:
java -jar MPC.jar hux_module.mpc
java -jar MPC.jar huxley_name_contains.mpc
java -jar MPC.jar Ado2Ino3comp.mpc
java -jar MPC.jar example.mpc
java -jar MPC.jar example2.mpc

2. Diff output with reference files:
 diff -s hux_module.mod ref/hux_module.mod
 diff -s huxley_name_contains.mod ref/huxley_name_contains.mod
 diff -s Ado2Ino3comp.mod ref/Ado2Ino3comp.mod
 diff -s example.mod ref/example.mod
 diff -s example2.mod ref/example2.mod

