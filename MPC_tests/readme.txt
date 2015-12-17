MPC tests

Generate JSim mod files and diff them with mod files in ref subdirectory 
to confirm no regressions in new code:

1. Generate JSim mod files:
java -jar MPC.jar hux_module.mpc
java -jar MPC.jar huxley_name_contains.mpc
java -jar MPC.jar Ado2Ino3comp.mpc
java -jar MPC.jar example.mpc
java -jar MPC.jar example2.mpc
java -jar MPC.jar test_WhiteSpace.mpc

2. Diff output with reference files:
diff -s ref/hux_module.mod hux_module.mod
diff -s ref/huxley_name_contains.mod huxley_name_contains.mod
diff -s ref/Ado2Ino3comp.mod Ado2Ino3comp.mod
diff -s ref/example.mod example.mod
diff -s ref/example2.mod example2.mod
diff -s ref/test_WhiteSpace.mod test_WhiteSpace.mod  
