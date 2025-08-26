# FRAMEWORK TESTS - testing the whole framework
macro(mirco_framework_test name_of_input_file)
  set (input_location ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file})
  add_test(NAME ${name_of_input_file}
  COMMAND ./mirco ${input_location})

endmacro(mirco_framework_test)

# List of tests
mirco_framework_test(input_sup2.yaml)
mirco_framework_test(input_sup5.yaml)
mirco_framework_test(input_sup6.yaml)
mirco_framework_test(input_sup7.yaml)
