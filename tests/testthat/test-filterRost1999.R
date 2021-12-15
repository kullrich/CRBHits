test_that("filter_rost1999() filters correct", {
    #artificial rbh
    rbh <- data.frame(rbind(
        c("x1", "y1", 5, 5),
        c("x1", "y1", 100, 5),
        c("x1", "y1", 15, 500),
        c("x1", "y1", 20, 500)))
    rbh.filtered <- filter_rost1999(rbh)
    expect_equal(dim(rbh.filtered)[1], 2)
})
