bookdown::render_book("index.Rmd", 
                      output_format = bookdown::gitbook(split_by = c("section+number")),
                      clean = TRUE, 
                      envir = parent.frame(), clean_envir = !interactive(), 
                      output_dir = "docs", 
                      new_session = NA, preview = FALSE, 
                      encoding = "UTF-8", config_file = "_bookdown.yml")