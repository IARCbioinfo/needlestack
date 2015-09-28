#!/usr/bin/awk -f
BEGIN {
  OFS = "\t"; # tab-delimited output
}
# Use substr instead of regex to match a starting ">"
substr($0, 1, 1) == ">" {
  if (seqlen) {
    # Only print info for this sequence if no target was given
    # or its id matches the target.
    if (! target || id == target) {
      print "##contig=<ID="id",length="seqlen">";
    }
  }
  # Get sequence id:
  # 1. Split header on whitespace (fields[1] is now ">id")
  split($0, fields);
  # 2. Get portion of first field after the starting ">"
  id = substr(fields[1], 2);
  seqlen = 0;
  next;
}
{
  seqlen = seqlen + length($0);
}
END {
  if (! target || id == target) {
    print "##contig=<ID="id",length="seqlen">";
  }
}
