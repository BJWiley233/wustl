### Just run this from the current directory with the command:

	R -e "shiny::runApp(launch.browser=T)"
	
	It should take about 10-20 seconds to load and read previous stored DESeq2 results.
	
	To use the app you just need to sart typing a gene name.  It will
	autocomplete and you can click the gene.  Sometimes it takes typing
	2 letters to get the gene names to populate.  To select a different
	gene, just click, delete, and type a new gene.
	
	There are a few more files here:
		1. First I cleaned the gene names because some were in Date form from Excel
		2. I ran the DESeq2 analysis in file 'task09.R' and wrote all results to
		   files so you don't have to wait for the analyis.  There is a cleaned
		   gene names with counts file 'counts_cleaned.txt' and three results files:
				a. Controlling for Gender     : 'controlled_gender_res.txt'
				b. Controlling for Age	      : 'controlled_age_res.txt'
				c. Controlling for Gender+Age : 'controlled_gender_age_res.txt'
		3. There are two tabs for Differential expression.  First tab shows
		   the DE for controlling factors individually testing Pathogenic Aging vs. Control.
		   The second tab shows DE for combining controlling terms and testing Pathogenic Aging vs. Control.
		  
