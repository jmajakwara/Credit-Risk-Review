###############################################################################
# Prosper P2P loans (Kaggle: henryokam/prosper-loan-data = prosperLoanData.csv, 113,937 rows x 81 variables)
###############################################################################
suppressPackageStartupMessages({ library(tidyverse) })

## ---- settings ---------------------------------------------------------------------------------
data_file  <- "DataProsper.csv"           # Kaggle file 
out_dir    <- "."
period_start <- as.Date("2009-07-31")         # keep loans originated AFTER this date (maintained)
default_statuses <- c("Defaulted", "Chargedoff", "Past Due (>120 days)")   # maintained
days_per_month   <- 30.4375
maturity_tol_days <- 31      # a completed loan closed >= (term - 31 days) after origination counts as paid to maturity
snapshot_override <- NULL    # e.g. as.Date("2014-03-12"); NULL = derive from the data (recommended)
bin_rule <- "mid"            # credit-score bin -> FICO band by bin "mid"point (reproduces the original mapping), "upper" or "lower"
prosperscore_11_as_na <- TRUE

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "prosper_cleaning_log.txt"), split = TRUE)
flow <- data.frame(step = character(), n = integer())
add_flow <- function(label, n) { flow[nrow(flow) + 1, ] <<- list(label, as.integer(n)); cat(sprintf("[%s] n = %d\n", label, n)) }

## ---- 1. read ---------------------------------------------------------------------------------------
raw <- read.csv(data_file, na.strings = c("NA", ""), stringsAsFactors = FALSE)     # names: ProsperRating..numeric. etc.
need <- c("LoanKey", "ListingKey", "LoanStatus", "ClosedDate", "LoanOriginationDate", "Term", "ProsperScore",
          "LoanFirstDefaultedCycleNumber", "LoanCurrentDaysDelinquent", "LoanMonthsSinceOrigination",
          "CreditScoreRangeLower", "CreditScoreRangeUpper", "FirstRecordedCreditLine")
stopifnot(all(need %in% names(raw)))
add_flow("rows in file", nrow(raw))

## ---- 2. duplicates (identical except ProsperScore) --------------------------------------------------------
n_key <- n_distinct(raw$LoanKey)
n_no_score <- nrow(distinct(select(raw, -ProsperScore)))
cat("distinct LoanKey:", n_key, "| distinct rows ignoring ProsperScore:", n_no_score,
    if (n_no_score == n_key) "-> duplicate copies differ ONLY in ProsperScore\n" else "-> WARNING: copies differ in other columns too\n")
conflict_keys <- raw %>% group_by(LoanKey) %>% summarise(k = n_distinct(ProsperScore, na.rm = TRUE) > 1, .groups = "drop") %>%
  filter(k) %>% pull(LoanKey)
cat("LoanKeys with conflicting ProsperScore among copies:", length(conflict_keys), "-> ProsperScore set to NA for them\n")
loans <- raw[!duplicated(raw$LoanKey), ]
loans$ProsperScore[loans$LoanKey %in% conflict_keys] <- NA
add_flow("unique loans (LoanKey)", nrow(loans))

## ---- 3. dates -----------------------------------------------------------------------------------------------
as_d <- function(x) as.Date(substr(x, 1, 10))
for (v in c("ListingCreationDate", "LoanOriginationDate", "ClosedDate", "FirstRecordedCreditLine", "DateCreditPulled"))
  loans[[v]] <- as_d(loans[[v]])

## ---- 4. data cut-off (snapshot) ---------------------------------------------------------------------------------
latest <- c(ClosedDate = max(loans$ClosedDate, na.rm = TRUE), ListingCreationDate = max(loans$ListingCreationDate, na.rm = TRUE),
            DateCreditPulled = max(loans$DateCreditPulled, na.rm = TRUE), LoanOriginationDate = max(loans$LoanOriginationDate, na.rm = TRUE))
cat("\nLatest dates in the file:\n"); print(as.Date(latest))
snapshot_date <- if (!is.null(snapshot_override)) snapshot_override else as.Date(max(latest))
act <- loans$LoanStatus %in% c("Current", "FinalPaymentInProgress")
impl <- loans$LoanOriginationDate + round(loans$LoanMonthsSinceOrigination * days_per_month)
cat("Cut-off used (snapshot_date):", format(snapshot_date), "\n")
cat("Cross-check: origination + LoanMonthsSinceOrigination for active loans has median", format(median(impl[act], na.rm = TRUE)),
    "(the old script used 2014-11-03, about 8 months too late)\n")

## ---- 5. period ---------------------------------------------------------------------------------------------------
loans <- loans %>% filter(LoanOriginationDate > period_start)
add_flow(paste("originated after", format(period_start)), nrow(loans))

## ---- 6. outcome, event type and time --------------------------------------------------------------------
known <- c("Current", "Completed", "Chargedoff", "Defaulted", "FinalPaymentInProgress", "Past Due (1-15 days)",
           "Past Due (16-30 days)", "Past Due (31-60 days)", "Past Due (61-90 days)", "Past Due (91-120 days)", "Past Due (>120 days)")
stopifnot(all(loans$LoanStatus %in% known), all(loans$Term %in% c(12, 36, 60)))
loans <- loans %>% mutate(
  term_days    = Term * days_per_month,
  elapsed_days = as.numeric(snapshot_date - LoanOriginationDate),
  closed_days  = as.numeric(ClosedDate - LoanOriginationDate),
  cycle_days   = LoanFirstDefaultedCycleNumber * days_per_month,
  is_default   = LoanStatus %in% default_statuses,
  is_completed = LoanStatus == "Completed",
  matured      = is_completed & !is.na(closed_days) & closed_days >= term_days - maturity_tol_days,
  event_type   = case_when(is_default ~ 1L, is_completed & !matured ~ 3L, TRUE ~ 2L),   # 1 default, 2 administrative, 3 early payoff
  time_days    = case_when(
    is_default & LoanStatus != "Past Due (>120 days)" & !is.na(cycle_days) ~ cycle_days,             # charge-off cycle
    is_default & LoanStatus != "Past Due (>120 days)"                      ~ closed_days,            # fallback: closing date
    is_default                                                             ~ elapsed_days - pmax(LoanCurrentDaysDelinquent - 120, 0),   # crossed 120 days past due
    event_type == 3L                                                       ~ closed_days,            # early payoff date
    is_completed                                                           ~ pmin(closed_days, term_days),   # paid to maturity
    TRUE                                                                   ~ pmin(elapsed_days, term_days)), # active at cut-off
  time_rule    = case_when(
    is_default & LoanStatus != "Past Due (>120 days)" & !is.na(cycle_days) ~ "charge-off cycle",
    is_default & LoanStatus != "Past Due (>120 days)"                      ~ "closed date",
    is_default                                                             ~ "120 days past due",
    event_type == 3L                                                       ~ "early payoff date",
    is_completed                                                           ~ "paid to maturity",
    TRUE                                                                   ~ "active at cut-off"))

cat("\nEvent type x status (1 default, 2 administrative, 3 early payoff):\n"); print(table(loans$LoanStatus, loans$event_type))
cat("\nTime rule used:\n"); print(table(loans$time_rule, loans$event_type))
d_both <- loans %>% filter(is_default, !is.na(cycle_days), !is.na(closed_days))
cat(sprintf("\nCharge-offs with both charge-off cycle and closed date: %d; |closed - cycle| > 60 days in %d (%.1f%%) -> cycle is used\n",
            nrow(d_both), sum(abs(d_both$closed_days - d_both$cycle_days) > 60), 100 * mean(abs(d_both$closed_days - d_both$cycle_days) > 60)))
bad <- is.na(loans$time_days) | loans$time_days <= 0
cat("Loans with missing or non-positive time (removed):", sum(bad), "\n"); print(table(loans$LoanStatus[bad], loans$event_type[bad]))
audit <- loans %>% select(LoanKey, LoanStatus, LoanOriginationDate, ClosedDate, Term, event_type, time_days, time_rule, elapsed_days, closed_days)
loans <- loans[!bad, ]
add_flow("positive follow-up time", nrow(loans))
loans <- loans %>% mutate(time = time_days,                                   # days 
                          status = factor(as.integer(event_type == 1L)))      # old 0/1 default indicator (0 now = both censoring types)

## ---- 7. derived design variables, drop identifiers / leakage / high missingness -------------------------------------
loans$OriginationDate  <- loans$LoanOriginationDate
loans$OriginationYear  <- factor(format(loans$LoanOriginationDate, "%Y"))
loans$CreditHistoryYears <- as.numeric(loans$LoanOriginationDate - loans$FirstRecordedCreditLine) / 365.25

id_var   <- c("ListingNumber", "ListingCreationDate", "ClosedDate", "GroupKey", "DateCreditPulled",
              "LoanFirstDefaultedCycleNumber", "LoanNumber", "MemberKey")                       # as before
leak_var <- c("LoanCurrentDaysDelinquent", "LoanMonthsSinceOrigination", grep("^LP_", names(loans), value = TRUE))   # outcome consequences
work_var <- c("term_days", "elapsed_days", "closed_days", "cycle_days", "is_default", "is_completed", "matured", "time_days", "time_rule")
high_miss <- c("TotalProsperLoans", "TotalProsperPaymentsBilled", "OnTimeProsperPayments", "ProsperPaymentsLessThanOneMonthLate",
               "ProsperPaymentsOneMonthPlusLate", "ProsperPrincipalBorrowed", "ProsperPrincipalOutstanding",
               "ScorexChangeAtTimeOfListing", "CreditGrade")                                    
rep_var  <- c("FirstRecordedCreditLine", "IncomeRange", "LoanOriginationDate", "LoanStatus",
              "LoanOriginationQuarter", "ProsperRating..Alpha.")                                 
miss_pct <- colMeans(is.na(loans)) * 100
cat("\nMissingness (%) of the variables dropped for high missingness, in the analysis period:\n"); print(round(miss_pct[high_miss], 1))
cat("Other variables above 50% missing in the period (should be only dropped/identifier ones):\n")
print(round(sort(miss_pct[miss_pct > 50], decreasing = TRUE), 1))
loans <- loans %>% select(-any_of(c(id_var, leak_var, work_var, high_miss, rep_var)))
cat("Removed outcome-leaking variables:", paste(leak_var, collapse = ", "), "\n")

## ---- 8. credit score: bin -> FICO band ---------------------------------------------------------------------------
# Prosper reports the score in 20-point bins [lower, lower+19].  Standard FICO bands:
#   Poor < 580 | Fair 580-669 | Good 670-739 | Very good 740-799 | Excellent >= 800
# Only the bin 660-679 straddles a boundary (669|670).  Its midpoint 669.5 lies above 669, so it is "Good" 
# Fair 600-659, Good 660-739, Very good 740-799, Excellent 800-899
stopifnot(!anyNA(loans$CreditScoreRangeLower), !anyNA(loans$CreditScoreRangeUpper),
          all(loans$CreditScoreRangeUpper - loans$CreditScoreRangeLower == 19))
ref_score <- switch(bin_rule, mid = (loans$CreditScoreRangeLower + loans$CreditScoreRangeUpper) / 2,
                    upper = loans$CreditScoreRangeUpper, lower = loans$CreditScoreRangeLower)
band_lv <- c("Poor", "Fair", "Good", "Very good", "Excellent")
loans$CreditScoreRange <- cut(ref_score, breaks = c(-Inf, 579, 669, 739, 799, Inf), labels = band_lv, ordered_result = TRUE)
stopifnot(!anyNA(loans$CreditScoreRange))
bins <- paste(loans$CreditScoreRangeLower, loans$CreditScoreRangeUpper, sep = "-")
orig <- dplyr::case_match(bins, c("600-619", "620-639", "640-659") ~ "Fair",
                          c("660-679", "680-699", "700-719", "720-739") ~ "Good",
                          c("740-759", "760-779", "780-799") ~ "Very good",
                          c("800-819", "820-839", "840-859", "860-879", "880-899") ~ "Excellent")
cat("\nCredit score: bins that the ORIGINAL mapping left unmatched (NA):", sum(is.na(orig)), "\n")
cat("New mapping identical to the original mapping for all matched bins:",
    all(as.character(loans$CreditScoreRange)[!is.na(orig)] == orig[!is.na(orig)]), "\n")
print(table(Bin = bins, Band = loans$CreditScoreRange))
loans$CreditScore <- (loans$CreditScoreRangeLower + loans$CreditScoreRangeUpper) / 2     # numeric midpoint (alternative to the bands)
loans$CreditScoreRange <- droplevels(loans$CreditScoreRange)                                # "Poor" does not occur (minimum 600)
cat("Levels kept:", paste(levels(loans$CreditScoreRange), collapse = " < "), "(Poor does not occur in the period)\n")
loans <- loans %>% select(-CreditScoreRangeLower, -CreditScoreRangeUpper)

## ---- 9. values that are really missing / out of range ---------------------------------------------------------------
n0 <- sum(loans$StatedMonthlyIncome == 0 & loans$EmploymentStatus %in% c("Employed", "Full-time", "Self-employed", "Part-time"), na.rm = TRUE)
loans$StatedMonthlyIncome <- ifelse(loans$StatedMonthlyIncome == 0 &
                                      loans$EmploymentStatus %in% c("Employed", "Full-time", "Self-employed", "Part-time"),
                                    NA, loans$StatedMonthlyIncome)       # zero income while working is impossible (zeros for 'Not employed' kept)
cat("\nStatedMonthlyIncome = 0 while employed set to NA:", n0, "\n")
cat("DebtToIncomeRatio == 0 set to NA:", sum(loans$DebtToIncomeRatio == 0, na.rm = TRUE), "| top-coded at 10.01 (flag DTI_topcoded):",
    sum(loans$DebtToIncomeRatio >= 10.01, na.rm = TRUE), "\n")
loans$DTI_topcoded <- loans$DebtToIncomeRatio >= 10.01
loans$DebtToIncomeRatio <- ifelse(loans$DebtToIncomeRatio == 0, NA, loans$DebtToIncomeRatio)
if (prosperscore_11_as_na) {
  cat("ProsperScore = 11 (documented scale is 1-10) set to NA:", sum(loans$ProsperScore == 11, na.rm = TRUE), "\n")
  loans$ProsperScore <- ifelse(loans$ProsperScore == 11, NA, loans$ProsperScore)
}

## ---- 10. names, types ------------------------------------------------------------------------------------------------------
loans <- loans %>% rename(ProsperRating = ProsperRating..numeric., ListingCategory = ListingCategory..numeric.,
                          TradesNeverDelinquent = TradesNeverDelinquent..percentage.)
loans$Term <- factor(loans$Term, levels = c(12, 36, 60))
loans$ProsperRating <- factor(loans$ProsperRating)
chr <- setdiff(names(loans)[sapply(loans, is.character)], c("LoanKey", "ListingKey"))
loans[chr] <- lapply(loans[chr], factor)                                   # IDs stay character
loans <- droplevels(loans)

first <- c("LoanKey", "ListingKey", "time", "status", "event_type", "BorrowerRate", "LenderYield", "ProsperScore",
           "EmploymentStatusDuration", "CurrentCreditLines", "OpenCreditLines", "TotalCreditLinespast7years", "OpenRevolvingAccounts",
           "OpenRevolvingMonthlyPayment", "InquiriesLast6Months", "TotalInquiries", "RevolvingCreditBalance", "CurrentDelinquencies",
           "AmountDelinquent", "DelinquenciesLast7Years", "PublicRecordsLast10Years", "PublicRecordsLast12Months", "BankcardUtilization",
           "AvailableBankcardCredit", "DebtToIncomeRatio", "LoanOriginalAmount", "PercentFunded", "StatedMonthlyIncome",
           "TradesNeverDelinquent", "IncomeVerifiable", "ProsperRating", "ListingCategory", "EmploymentStatus", "Occupation",
           "CurrentlyInGroup", "IsBorrowerHomeowner", "CreditScoreRange", "CreditScore", "Term", "OriginationDate", "OriginationYear",
           "CreditHistoryYears")
loans <- loans %>% select(any_of(first), everything())
add_flow("final analysis file (all terms)", nrow(loans))

## ---- 11. summaries ---------------------------------------------------------------------------------------------------------------
cat("\n== Event types by term (1 default, 2 administrative, 3 early payoff) ==\n"); print(table(loans$Term, loans$event_type))
cat("\n== Event types by origination year ==\n"); print(table(loans$OriginationYear, loans$event_type))
cat("\n== Follow-up time (days) by event type ==\n"); print(tapply(loans$time, loans$event_type, summary))
cat("\n== Missingness (%) in the analysis file, variables with any missing ==\n")
mp <- colMeans(is.na(loans)) * 100; print(round(sort(mp[mp > 0], decreasing = TRUE), 2))

## ---- 12. write ---------------------------------------------------------------------------------------------------------------------------

write.csv(loans, file="prosper_analysis.csv", row.names = FALSE)
