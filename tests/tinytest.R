if ( requireNamespace("tinytest", quietly=TRUE) ){
  options(x12.delete = TRUE)
  tinytest::test_package("simPop")
}
