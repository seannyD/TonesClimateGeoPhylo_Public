library(dplyr)
library(tidyr)
library(purrr)

bt_read.log <- function(filename){
  con = file(filename, "r")
  i = 1
  j = 1
  attr_list = list()
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if (grepl("Iteration\t", line)|grepl("Tree\t", line)||grepl("Tree No\t", line)) {
      break
    }
    i = i + 1
  }
  close(con)
  
  settings_raw = readLines(filename, n = i - 1)
  settings = .get_attributes(settings_raw)
  
  d = read.table(filename, skip = i-1, sep = '\t',
                 header = TRUE, check.names = FALSE, quote=NULL)
  d[sapply(d, function(x) all(is.na(x)))] <- NULL
  
  attr(d,"settings") <- settings
  class(d) <- append(class(d), c("bt_log"))
  d
}

.get_attributes = function(line){
  one = dplyr::as_tibble(line)
  two = tidyr::separate(one, value, into = c("header", "info"), sep = "\\s{2,}", extra = "merge", fill = "right")
  three = dplyr::mutate(two, header = dplyr::na_if(header, ""))
  four = tidyr::fill(three, header) # fills empty info
  five = dplyr::filter(four, !is.na(info)) # gets rid of titles with empty info
  six = dplyr::mutate(five, info = stringr::str_trim(info))
  seven = split(six, six$header)
  purrr::map(seven, dplyr::pull, info)
}