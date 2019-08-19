#! R

##read all sheets into list of tibbles
library(readxl)
read_sheets_to_list <- function(FILENAME, COLTYPES = NULL) {

    ##allow coltypes input
    if(is.null(COLTYPES)){
      COLTYPES <- c()
    }

    ##define alll sheets, could be modified to match/grep on an input variable
    sheets <- readxl::excel_sheets(FILENAME)

    ##don't want empty sheets which excel helpfully initiates in workbooks
    sheetsLength <- unlist(lapply(sheets, function(f){
       sf <- readxl::read_excel(FILENAME, sheet = f)
       dim(sf)[1]
    }))

    ##for coltype date input need to find colnames of DOB, Date*
    sheetsDate <- lapply(sheets, function(f){
       sf <- readxl::read_excel(FILENAME, sheet = f)
       coltypes <- ifelse(colnames(sf) %in% c("DOB", grep("Date", colnames(sf), value=TRUE)), "date", "guess")
    })
    names(sheetsDate) <- sheets

    ##return the list of non-empty sheets, named as per original
    simpledate <- function(d){d <- as.Date(d); format(d, "%d-%m-%Y")}
    sheetsList <- lapply(sheets[sheetsLength>0], function(f){
       readxl::read_excel(FILENAME, sheet = f, col_types = sheetsDate[[f]]) %>%
       dplyr::mutate_at(vars(grep("date", sheetsDate[[f]])), simpledate)
    })
    names(sheetsList) <- sheets[sheetsLength>0]
    return(sheetsList)
}
