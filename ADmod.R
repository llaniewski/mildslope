#!/usr/bin/env Rscript

library(optparse)
options <- list(
        make_option(c("-f","--file"), "store", default="", help="Input file", type="character"),
        make_option(c("-o","--out"), "store", default="", help="Output file", type="character"),
	make_option(c("-x","--fix"), "store", default="", help="variables to fix", type="character")
)

opt <- parse_args(OptionParser(usage="Usage: ADmod -f inputfile [-o outputfile]", options))

if (opt$file == "") stop("Input file not specified\nUsage: ADmod -f file\n");
if (opt$out == "") { opt$out = paste(opt$file, "_",sep="") }

fix = strsplit(opt$fix," ")[[1]]
cat("To fix:\n")
print(fix);


f = file(opt$file)

lines = readLines(f)
close(f)


bracel = gregexpr("[{]",lines)
bracer = gregexpr("[}]",lines)

a = sapply(bracel,function(x){sum(x>0)})
b = sapply(bracer,function(x){sum(x>0)})
a = cumsum(a-b)

a[a>1]=1
begins = which(diff(a)==1)+1

f = file(opt$out)
open(f,"wt")
pushi = grep("push([rR]eal|[cC]ontrol)",lines)
looki = grep("look([rR]eal|[cC]ontrol)",lines)
popi = grep("pop([rR]eal|[cC]ontrol)",lines)

begins = c(begins,length(lines))
alli = sort(c(pushi,popi,begins,looki))
idx = 0
tmpname = "keep";
si = 0
buf = c()
vars = c()
decl = 0
control_stack = FALSE
for (i in alli) {
	if (i-si > 1) {
		buf = c(buf, lines[(si+1):(i-1)] );
	}
	if (i %in% begins) {
		buf = c(buf, lines[i] );
#		cat(vars,sep="\n")
#		cat(buf,sep="\n")
		if (length(vars) > 0) writeLines(text=vars,con=f)
		writeLines(text=buf,con=f)

		vars = c()
		buf = c()
		decl = 0;
	} else {
		l = lines[i]
		res = regmatches(l, regexec(
			"([[:space:]]*)((pop|push|look)([Rr]eal|[Cc]ontrol)([0-9][0-9]*)([^(]*)\\(([&]?)([^)]*)\\))",
			l))[[1]]
		
		tp = NULL
		ar = FALSE
		if (res[5] == "Real") {
			if (res[6] == "4") tp = "float"
			if (res[6] == "8") tp = "double"
		} else if (res[5] == "Control") {
			tp = "int";
		}
		if (res[7] == "array") ar = TRUE
		if (is.null(tp)) stop("Unknown type of push/pop: ",l);
		l1 = res[9]
		if (ar) {
			l2 = sub("^[^,]*,[ ]*","", l1);
			ar_size = as.integer(l2);
			ar_dec = paste("[",ar_size,"]",sep="")
			ar_idx = paste("[",1:ar_size-1,"]",sep="")
			l1 = sub(",.*$","",l1);
		} else {
			ar_dec = ""
			ar_idx = ""
		}
		var = l1;
		space = res[2]
		com = paste("// ADmod.R: ", res[3], sep="")
		if (var %in% fix) {
			cat("var: ",var," ----- fixed\n");
			buf = c(buf, com);
		} else if (res[5] == "Control") {
			name = "control_stack"
			bits = as.integer(res[6])
			mask = sprintf("0x%04x",2^bits-1)
			if (grepl("push", l)) {
				if (!control_stack) {
					vars = c(vars, paste(tp," ",name, " = 0x0000; ",com,sep=""));
					control_stack = TRUE
				}
				buf = c(buf, paste(space,name," = (", name, "<<", bits, ") + ", var, "; ",com,sep=""));
			} else if (grepl("look", l)) {
				buf = c(buf, paste(space,var, " = ", name, " & ",mask,"; ",com,sep=""));
			} else if (grepl("pop", l)) {
				buf = c(buf, paste(space, "{ ", var, " = ", name, " & ",mask,"; ",
					name, " = ", name, " >> ", bits,"; } ",com,sep=""));
			} else {
				stop("Unknown type of push/pop: ",l);
			}
		} else {
			if (grepl("push", l)) {
				idx = idx + 1
				name = paste(tmpname, idx, sep="_")
				if (idx > decl) {
					vars = c(vars, paste(tp," ",name,ar_dec,"; ",com,sep=""));
					decl = idx;
				}
				buf = c(buf, paste(space,name,ar_idx," = ", var,ar_idx, "; ",com,sep=""));
			} else if (grepl("look", l)) {
				name = paste(tmpname, idx, sep="_")
				buf = c(buf, paste(space,var,ar_idx, " = ", name,ar_idx,"; ",com,sep=""));
			} else if (grepl("pop", l)) {
				name = paste(tmpname, idx, sep="_")
				buf = c(buf, paste(space,var,ar_idx, " = ", name,ar_idx,"; ",com,sep=""));
				idx = idx - 1;
			} else {
				stop("Unknown type of push/pop: ",l);
			}
		}
	}
	si = i;
}
close(f)
