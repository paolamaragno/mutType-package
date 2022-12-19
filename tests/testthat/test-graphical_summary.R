library('pdftools')

test_that("Graphical visualization of the mutation types with a frequency higher than a threshold", {

  c_table <- data.frame(mut_types = c("AA[C>A]AA", "AA[C>A]AC", "AA[C>A]AG", "AA[C>A]AT", "AA[C>A]CA", "AA[C>A]CC",
                                     "AA[C>A]CT", "AA[C>A]GA", "AA[C>A]GC", "AA[C>A]GG", "AA[C>A]GT", "AA[C>A]TA",
                                     "AA[C>A]TC"),
                        Freq = c(15,  13,  32,  44,  25,  5,  67,  21,  14,  42,  52,  21,  19))

  c <- c_table[c_table$Freq >= 30,]

  x <- c$mut_types

  f <- gl(ceiling(nrow(c)/100), 100, length=nrow(c))

  pdf("Mut_types_visualization_30_test.pdf")

  tapply(x, f, FUN = function(sublist){c2 <- c[which(c$mut_types %in% sublist),]
        print(ggplot(c2, aes(x = mut_types, y = Freq)) +
                geom_bar(stat = "identity") +
                geom_col(aes(fill = Freq)) +
                scale_fill_gradient2(low = 'yellow',
                                     mid ='red',
                                     high = 'darkgreen',
                                     midpoint = max(c2$Freq)/2) +
                theme_classic() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 5)) +
                labs(title = paste("Frequency of each mutation type with Freq >", as.character(30)),
                     subtitle = paste('Overall number of mutation types with Freq >', as.character(30), ':',  as.character(nrow(c)),
                     '\nOverall number of mutation types with Freq >', as.character(30), 'in this page:',  as.character(nrow(c2))),
                     x = "Mutation type",
                     y = "Frequency"))
  })

  dev.off()

  expect_equal(pdf_data(graphical_summary(c_table, 30, "Mut_types_visualization_30.pdf")), pdf_data("Mut_types_visualization_30_test.pdf"))
})
