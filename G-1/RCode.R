#
# Thermoanaerobacter kivui G1 genome curation and annotation
#
# Author: GG. Thallinger, Graz University of Technology, 2024


# =============================================================================
# Start of loading packages


# =============================================================================
# Start of functions

# Write a progress message to the console and flush it for immediate display
#
writeLog <- function(textstr, verbose=TRUE, continue=FALSE)
{
  if ((class(textstr) == "character") && (length(textstr) == 1)) {
    if (continue) { # does not advance to the next line in the output
      output_string <- textstr
    } else {
      output_string <- paste0(textstr, "\n")
    }
    if (verbose) {
      cat(output_string)
    }
  } else {
    if (verbose) {
      output_string <- print(textstr)
    } else {
      output_string <- ""
    }
  }
  if (verbose) {
    flush.console()
  }
  return(output_string)
}


# Return start/end line indices and feature name of the feature
# the current line index belongs to
#
getFeatureStartEnd <- function(lt.idx, gbk.text)
{
  curr.idx <- lt.idx
  repeat { # Search upwards
    if ((length(grep("^     [A-Za-z_]+ ", gbk.text[curr.idx])) == 1) | # Found feature
        (length(grep("^[A-Z]+", gbk.text[curr.idx])) == 1)) {          # Found section, emergency break
      break  # New feature or section
    }
    curr.idx <- curr.idx - 1
  }
  start.idx <- curr.idx
  curr.idx <- lt.idx
  repeat { # Search downwards
    curr.idx <- curr.idx + 1
    if ((length(grep("^     [A-Za-z_]+ ", gbk.text[curr.idx])) == 1) | # Found feature
        (length(grep("^[A-Z]+", gbk.text[curr.idx])) == 1)) {          # Found section, emergency break
      break  # New feature or section
    }
  }
  end.idx <- curr.idx - 1
  feature <- sub("^     ([A-Za-z_]+) .*", "\\1", gbk.text[start.idx])
  return(list(Start=start.idx, End=end.idx, Feature=feature))
}


# Return the line index of a specific GenBank locus_tag, feature and 
# tag combination; if not present return 0
#
getTagIndex <- function(locus.tag, feature, tag, gbk.text)
{
  feature.idx <- 0
  feature.name <- ""
  tag.idx <- 0
  lt.idx <- grep(paste0('    /locus_tag="', locus.tag, '"'), gbk.text)
  if (length(lt.idx) == 0) {
    writeLog(sprintf("getTagIndex: /locus_tag '%s' not found", locus.tag))
  } else {
    features <- do.call(rbind.data.frame, lapply(lt.idx, getFeatureStartEnd, gbk.text))
    ft.idx <- grep(feature, features$Feature)
    if (length(ft.idx) == 0) { # Specific feature not found
      writeLog(sprintf("getTagIndex: Feature '%s' for /locus_tag '%s' not found", feature, locus.tag))
    } else {
      lt.idx <- lt.idx[ft.idx]
      feature.idx <- features$Start[ft.idx]
      tag.idx <- grep(tag, gbk.text[feature.idx:features$End[ft.idx]])
      if (length(tag.idx) == 0) { # Specific feature not found
        tag.idx <- 0
      } else {
        tag.idx <- tag.idx + features$Start[ft.idx] - 1
        feature.name <- features$Feature[ft.idx]
      }
    }
  }
  return(list(LTidx=lt.idx, Fidx=feature.idx, Tidx=tag.idx))
}


# Return all line indices belonging to a tag (start and continuation lines)
# 
getTagLines <- function(curr.idx, gbk.text)
{
   cont.idx <- NULL
   repeat {
     cont.idx <- c(cont.idx, curr.idx)
     curr.idx <- curr.idx + 1
     if ((length(grep("^(  |FT)                   /", gbk.text[curr.idx])) == 1) | 
         (length(grep("^(  |FT)   [A-Za-z_]+ ", gbk.text[curr.idx])) == 1) | 
         (length(grep("^SQ   ", gbk.text[curr.idx])) == 1) | 
         (length(grep("^ORIGIN", gbk.text[curr.idx])) == 1)) { 
       break  # New tag, feature or section
     }
   }
   return(cont.idx)
}


# Return the line indices of a specific GenBank tag and 
# all subsequent lines beloning to the tag
#
getTagIndices <- function(tag, gbk.text)
{
  tag.idx <- grep(tag, gbk.text)
  writeLog(sprintf("getTagIndices: %.40s   %d", tag, length(tag.idx)))
  all.idx <- unlist(sapply(tag.idx, getTagLines, gbk.text))
  writeLog(sprintf("getTagIndices: %.40s   %d", tag, length(all.idx)))
  return(all.idx)
}


# Paste GenBank tags to the given string and wrap lines if necessary
# Return the strings in the specified number of lines
#
pasteTag <- function(locus.tag, prefix, tag, value, olines, psplit="(?<=.{58})")
{
  tag.full <- sprintf('%s%s', tag, value)
  if (nchar(tag.full) <= 58) {
    tag.split <- tag.full
  } else {
    if (psplit == " ") {
      sstart <- 1
      send <- NULL
      split.pos <- unlist(gregexpr(" ", paste0(tag.full, " ")))
      if ((length(split.pos) == 0) | (split.pos[1] > 58)) {
        warning(sprintf("No valid split positions: %s", paste0(split.pos, collapse=",")))
        tag.split <- unlist(strsplit(tag.full, "(?<=.{58})" , perl=TRUE))
      } else {
        split.lim <- 58
        for (i in seq_len(length(split.pos))) {
          if (split.pos[i] > split.lim) {
            send <- c(send, split.pos[i-1] - 1)
            sstart <- c(sstart, split.pos[i-1] + 1)
            split.lim <- split.pos[i-1] + 1 + 58
          }
        }
        send <- c(send, nchar(tag.full))
        tag.split <- substring(tag.full, sstart, send)
      }
    } else {
      tag.split <- unlist(strsplit(tag.full, psplit, perl=TRUE))
    }
  }
  if (length(tag.split) < olines) {
    warning(sprintf("%s - Requested number of output lines (%d) smaller than number of tag lines (%d)",
            locus.tag, olines, length(tag.split)))
    tag.split <- c(tag.split, rep(" ", olines - length(tag.split)))
  }
  tag.lines <- sprintf("%21.21s%s", " ", tag.split) # Add spacing before tag lines
  tag.out <- tag.lines[1:olines]
  if (olines < length(tag.lines)) { # Paste the remaining lines to the last output line
    tag.add <- paste0("\n", tag.lines[(olines+1):length(tag.lines)])
    tag.out[olines] <- paste0(tag.out[olines], do.call(paste0, as.list(tag.add)))
  }
  if (prefix != "") {
    tag.out[1] <- paste0(prefix, "\n", tag.out[1])
  }
  return(tag.out)
}


# Update the annotation (in GenBank format) derived by transfer of an 
# existing one by Genious and integrate manual annotations present in
# an tab separated file:
#   - Remove all comments from Genious
#   - Update /gene and /product texts using contents of the file
#   - Change /locus_tags to the specified prefix and number them 
#     continuously with the given offset.
#
#  Genious introduced tags:
#    /note="Derived using Geneious Prime 2023.2.1 'Annotate
#    from Database' based on nucleotide similarity"
#    /Transferred_From="NZ_CP009170"
#    /Transferred_Similarity="100.00%"
#    /Primary_Match="100.00%; Adjusted-Interval-0:0/0;
#    carbohydrate ABC transporter permease CDS (<a
#    href=""urn:local:.:i-h9q53gy"">NZ_CP009170</a>)"
#    /label="carbohydrate ABC transporter permease CDS"
#    /Primary_Match="100.00%; Adjusted-Interval-0:0/0; repeat
#    region (<a href=""urn:local:.:i-h9q53gy"">NZ_CP009170</a>)"
#
updateGenBankAnnotation <- function(gbk.name, uname, prefix="TKV_G1", loffset=5, out.name, verbose=FALSE)
{
  gbk <- readLines(gbk.name)

  # Some statistics
  ts.idx <- grep('/Transferred_Similarity=', gbk)
  ts.100.idx <- grep('/Transferred_Similarity="100.00%', gbk)
  ts.diff.idx <- setdiff(ts.idx, ts.100.idx)
  writeLog(sprintf("Transferred: %5d   identical: %5d    differing: %5d",
                   length(ts.idx), length(ts.100.idx), length(ts.diff.idx)))

  # Remove tags and changes introduced by Genious
  to.remove <- c('/note="Derived using Geneious Prime', '/Transferred_From="',
                 '/Transferred_Similarity="', '/Primary_Match="', '/label="',
                 '/Alternative_Match="', 'Adjusted_Transferred_Similarity="', '/old_locus_tag="')
  rem.idx <- unlist(sapply(to.remove, getTagIndices, gbk))
  rem.idx.unique <- unique(rem.idx)
  gbk.red <- gbk[-rem.idx]

  gbk.red <- sub("/Transferred_Translation=", "/translation=", gbk.red)

  # Update the annotation using the supplied file, extract loci 
  # with annotation to be updatad and remove duplicate locus_tags
  ann.upd <- read.table(uname, header=TRUE, sep="\t", na="")
  writeLog(sprintf("Annotations: %5d    Short: %5d    Long %5d",
                   nrow(ann.upd), length(which(!is.na(ann.upd$gene))),
                   length(which(!is.na(ann.upd$product)))))
  for (i in seq_len(nrow(ann.upd))) {
    lt <- trimws(ann.upd$locus_tag[i])
    gene <- trimws(ann.upd$gene[i])
    product <- trimws(ann.upd$product[i])
    if (!is.na(gene)) { # update /gene tag
      # Update gene tag for "gene" feature
      tag.info <- getTagIndex(lt, "gene", "/gene", gbk.red)
      if (tag.info$Fidx > 0) {     # Locus tag and feature were found
        if (tag.info$Tidx > 0) {   # /gene tag is present, replace content
          old.gene <- sub('.*/gene="([A-Za-z0-9]+)"', "\\1", gbk.red[tag.info$Tidx])
          gbk.red[tag.info$Tidx] <- sub(old.gene, gene, gbk.red[tag.info$Tidx])
          writeLog(sprintf("%s - Replaced short name for '%s' from: '%s' to '%s'", lt, "gene", old.gene, gene), verbose=verbose)
        } else {
          gbk.red[tag.info$Fidx] <- paste0(gbk.red[tag.info$Fidx], '\n                     /gene="', gene, '"')
          writeLog(sprintf("%s - Added short name for '%s': '%s'", lt, "gene", gene), verbose=verbose)
        }
      }
      # Update gene tag for "CDS" feature
      tag.info <- getTagIndex(lt, "CDS", "/gene", gbk.red)
      if (tag.info$Fidx > 0) {     # Locus tag and feature were found
        if (tag.info$Tidx > 0) {   # /gene tag is present, replace content
          old.gene <- sub('.*/gene="([A-Za-z0-9]+)"', "\\1", gbk.red[tag.info$Tidx])
          gbk.red[tag.info$Tidx] <- sub(old.gene, gene, gbk.red[tag.info$Tidx])
          writeLog(sprintf("%s - Replaced short name for '%s' from: '%s' to '%s'", lt, "CDS", old.gene, gene), verbose=verbose)
        } else {
          gbk.red[tag.info$Fidx] <- paste0(gbk.red[tag.info$Fidx], '\n                     /gene="', gene, '"')
          writeLog(sprintf("%s - Added short name for '%s': '%s'", lt, "CDS", gene), verbose=verbose)
        }
      }
    }
    if (!is.na(product)) { # update /product tag
      tag.info <- getTagIndex(lt, "CDS", "/product", gbk.red)
      if (tag.info$Fidx > 0) {     # Locus tag and feature were found
        if (tag.info$Tidx > 0) {   # /product tag is present, replace content
          old.product <- sub('.*/product="([^"]+).*', "\\1", gbk.red[tag.info$Tidx])
          prod.lines <- getTagLines(tag.info$Tidx, gbk.red)
          writeLog(sprintf("%s - Replaced product from: '%s' [%d] to '%s'", lt, old.product, 
                           length(prod.lines), product), verbose=verbose)
          pasted.lines <- pasteTag(lt, "", '/product="', paste0(product, '"'), length(prod.lines), " ")
          gbk.red[prod.lines] <- pasted.lines
          if (lt == "TKV_RS00355x") {
            print(prod.lines)
            print(pasted.lines)
            print(gbk.red[prod.lines])
          }
        } else {
          writeLog(sprintf("%s - Added product name: '%s'", lt, product), verbose=verbose)
          stop(sprintf("%s - Adding of /product not supported", lt))
        }
      }
    }
  }

  # Add a chromosome tag (may be checked by the validator)
  orga.idx <- grep("/organism=", gbk.red)
  gbk.red[orga.idx] <- paste0(gbk.red[orga.idx], '\n                     /chromosome="1"')

  # Remove any empty lines introduced in any of the previous steps
  empty.idx <- grep(sprintf("^%21.21s $", ""), gbk.red)
  if (length(empty.idx)) {
    gbk.red <- gbk.red[-empty.idx]
  }

  # Finally change the /locus_tag to /old_locus_tag and add a new consecutivly 
  # numbered /locus_tag
  lt.idx <- grep(' /locus_tag=', gbk.red)
  lt.table <- table(gbk.red[lt.idx])
  writeLog(sprintf("/locus_tag entries: %5d   Pairs: %5d    Multiples: %5d     Other: %5d",
                   length(lt.idx), length(lt.idx/2), length(which(lt.table == 2)),
                   length(which(lt.table > 2))))

  gbk.red[lt.idx] <- sub(" /locus_tag=", " /old_locus_tag=", gbk.red[lt.idx])
  if (length(lt.idx) %% 2 != 0) {
    stop(sprintf("Number of /locus_tag entries not a multiple of 2"))
  }
  # FIXXXME: assumes that tags for the same locus are in pairs and consecutive.
  lt.defined <- length(lt.idx) / 2
  rnum <- rep(sprintf("%05d", 1:lt.defined * loffset), each=2)
  gbk.red[lt.idx] <- paste0('                     /locus_tag="', prefix, rnum, '"\n', gbk.red[lt.idx])

  print(c(sub(".gb$", "", basename(gbk.name)), sub(".gb$", "", basename(out.name))))
  acc.old <- sub("LOCUS +([^ ]+) .*", "\\1", gbk.red[1])
  acc.new <- sub(".gb$", "", basename(out.name))
  gbk.red[1:3] <- sub(acc.old, acc.new, gbk.red[1:3])
  cdate <- toupper(format(Sys.time(), "%e-%b-%Y"))
  gbk.red[1] <- sub("[0-9]+-[A-Z]+-[0-9]+$", cdate, gbk.red[1])
  print(cat(gbk.red[1:3], sep="\n"))

  adate.txt <- "            Annotation Date                   :: "
  adate.idx <- grep(adate.txt, gbk.red)
  gbk.red[adate.idx] <-  paste0(adate.txt, toupper(format(Sys.time(), "%e-%b-%Y %H:%M:%S")))

  writeLines(gbk.red, out.name)
  writeLog(sprintf("Initial length: %6d   - Reduced length: %6d",
                   length(gbk), length(gbk.red)))
  return(gbk.red)
}


# =============================================================================
# Start of analysis
start_analysis <- function() {}

gname <- file.path(adir, "G1_Canu_polished.gb")
uname <- file.path(adir, "G1_Annotation-Updates.tsv")
oname <- file.path(adir, "G1_Canu_updated.gb")
g1.upd <- updateGenBankAnnotation(gname, uname, prefix="TKVG1_", out.name=oname)

# End of analysis
end_analysis <- function() {}
