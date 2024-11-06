test_that("filter_rost1999() filters correct", {
    #artificial rbh
    rbh <- data.frame(
        query_id = rep("x1", 4),
        subject_id = rep("y1", 4),
        perc_identity = c(5, 100, 15, 20),
        alignment_length = c(5, 5, 500, 500)
    )
    rbh.filtered <- filter_rost1999(rbh)
    expect_equal(dim(rbh.filtered)[1], 2)
})
