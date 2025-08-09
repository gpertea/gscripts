#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

#-------------------------------------------------------------------------------
#   Usage
#-------------------------------------------------------------------------------
#
#   This script removes multi-line C-style (/* ... */) comments from a
#   C, C++, or Java source file if they span a given number of lines.
#
#   SYNTAX:
#   perl remove_comments.pl [-l N] < input_file > output_file
#
#   ARGUMENTS:
#   -l N : Optional. The minimum number of lines a comment block must
#          span to be removed. Defaults to 4.
#
#   BEHAVIOR:
#   - Reads from STDIN and prints to STDOUT.
#   - Correctly ignores '/*' or '*/' if they appear after a '//' single-line
#     comment marker on the same line.
#   - If a multi-line comment starts after some code on a line, that code
#     is preserved.
#   - If a multi-line comment ends before some code on a line, that code
#     is also preserved.
#   - Unclosed comments are preserved at the end of the file to prevent
#     accidental code loss.
#
#-------------------------------------------------------------------------------

# --- Configuration ---
my %opts;
getopts('l:', \%opts);
my $MIN_LINES = $opts{l} || 4;

# --- State Variables ---
my $in_multiline_comment = 0;   # Flag: 1 if we are inside a /* ... */ block, 0 otherwise.
my @comment_buffer;             # Stores the lines of the current comment block being evaluated.
my $comment_start_prefix = "";  # Stores any code that appeared *before* the '/*' on the starting line.

# --- Main Loop ---
# Process the input file (from STDIN) line by line.
while (my $line = <>) {
    if ($in_multiline_comment) {
        # STATE: We are currently inside a multi-line comment block.
        push @comment_buffer, $line;

        # Check if this line closes the comment block.
        # Note: A '*/' cannot be commented out by '//' on the same line
        # because the '/*' on a previous line has already put us into this state.
        if ($line =~ /\*\//) {
            # Found the end of the comment block: '*/'
            my $code_after_comment = $'; # The part of the string AFTER the '*/' match.

            if (scalar(@comment_buffer) >= $MIN_LINES) {
                # The block is long enough to be removed.
                # We stitch the prefix from the start line with the postfix from this end line.
                my $stitched_line = $comment_start_prefix . $code_after_comment;

                # Print the stitched line only if it contains non-whitespace characters.
                # This completely removes comment blocks that occupied whole lines,
                # rather than replacing them with a blank line.
                if ($stitched_line =~ /\S/) {
                    print $stitched_line;
                }
            } else {
                # The block is too short, so we print the buffered content unchanged.
                print @comment_buffer;
            }

            # Reset state for the next lines.
            $in_multiline_comment = 0;
            @comment_buffer = ();
            $comment_start_prefix = "";
        }
    } else {
        # STATE: We are not in a multi-line comment; we are looking for a start.
        # Check for a potential start '/*'
        if ($line =~ /\/\*/) {
            my $potential_start_pos = $-[0]; # Get the start position of the '/*' match.
            my $is_real_start = 1;

            # Now, validate that this '/*' is not itself commented out by '//'.
            # We check if a '//' marker exists on the line *before* our '/*' match.
            if ($line =~ m|//|) {
                my $single_line_comment_pos = $-[0]; # Get start position of '//' match.
                if ($single_line_comment_pos < $potential_start_pos) {
                    # The '//' comes first, so the '/*' is part of a single-line comment.
                    # This is a false start, so we should ignore it.
                    $is_real_start = 0;
                }
            }

            if ($is_real_start) {
                # This appears to be a real start of a multi-line comment.
                my $prefix = substr($line, 0, $potential_start_pos);
                my $rest_of_line = substr($line, $potential_start_pos + 2); # +2 to skip the '/*' itself.

                # Check if the comment also ENDS on this same line.
                if ($rest_of_line =~ /\*\//) {
                    # It's a single-line block comment, e.g., int x = 1; /* comment */
                    # We do nothing to it, just print the original line and move on.
                    print $line;
                } else {
                    # This is the start of a true multi-line comment.
                    # Change state and start buffering.
                    $in_multiline_comment = 1;
                    push @comment_buffer, $line;
                    $comment_start_prefix = $prefix; # Save any code that came before '/*'.
                }
            } else {
                # The '/*' was found, but it was inside a '//' comment. Just print the line.
                print $line;
            }
        } else {
            # This line does not start a multi-line comment. Just print it.
            print $line;
        }
    }
}

# --- End of File ---
# If the file ends while we are still inside a comment block, it's an unclosed
# comment. We should print the buffer to avoid losing code or breaking the file.
if ($in_multiline_comment) {
    print @comment_buffer;
}
