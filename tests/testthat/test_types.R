context("DGEobj - tests for types.R functions")


test_that('types.R: baseType()', {
    expect_equal(baseType(DGEobj, "intensity"), "assay")
    expect_equal(baseType(DGEobj, "design"), "col")
    expect_equal(baseType(DGEobj, "intensity_orig"), "meta")

    expect_equivalent(baseType(DGEobj, "counts"), "assay")

    expect_error(baseType(DGEobj, "dog"))
})

test_that('types.R: baseTypes()', {
    expect_setequal(baseTypes(), c("row", "col", "assay", "meta" ))
    expect_setequal(baseTypes(DGEobj), c("row", "col", "assay", "meta"))
})

test_that('types.R: showTypes()', {
    showTypes_DGEobj <- showTypes(DGEobj)

    expect_s3_class(showTypes_DGEobj, "knitr_kable")
    expect_equal(length(showTypes_DGEobj), 50)

    warning('function does not appear to be working')
    # showTypes_notPretty_DGEobj <- showTypes(DGEobj, pretty = FALSE)
    # expect_equal(length(showTypes_notPretty_DGEobj), 50)
})

test_that('types.R: newType()', {
    newType_DGEobj <- newType(DGEobj, "MyType", "meta")
    expect_true("MyType" %in% names(attr(newType_DGEobj, "objDef")$type))
    expect_equal(attr(newType_DGEobj, "objDef")$type[["MyType"]], "meta")
    expect_false("MyType" %in% attr(newType_DGEobj, "objDef")$uniqueType)

    newType_DGEobj <- newType(DGEobj, "MyType", "assay", uniqueItem = TRUE)
    expect_true("MyType" %in% names(attr(newType_DGEobj, "objDef")$type))
    expect_equal(attr(newType_DGEobj, "objDef")$type[["MyType"]], "assay")
    expect_true("MyType" %in% attr(newType_DGEobj, "objDef")$uniqueType)
})

test_that('types.R: incorrect usage', {
    expect_error(DGEobj::showTypes())

    expect_error(DGEobj::newType())
    expect_error(DGEobj::newType(DGEobj))
    expect_error(DGEobj::newType(DGEobj, "MyType", "badType"))
})
