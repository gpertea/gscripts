#!/usr/bin/perl
use LWP::UserAgent;
use HTML::TreeBuilder;
use HTML::FormatMarkdown;

my $url = 'https://www.cog-genomics.org/plink/2.0/';
# get the base path after the domain name
my ($basepath) = $url =~ m{^https?://[^/]+(/.*)$};
my $ua = LWP::UserAgent->new;
my $response = $ua->get($url);

die "Error: ".$response->status_line."\n" unless $response->is_success;

my $tree = HTML::TreeBuilder->new_from_content($response->decoded_content);
# Find the td element with id='leftbar_full' and then extract links within it
my $leftbar = $tree->look_down(_tag => 'td', id => 'leftbar_full');
die("Error: couldn't find leftbar_full\n") unless $leftbar;
my $markdown_content = "";
append_md(\$markdown_content, $response); # Append to markdown content
# Find all unique links within the main page to scrape
my %unique_links;

my @links = $leftbar->look_down(_tag => 'a', href => qr{^\Q$basepath\E}); 
foreach my $link (@links) {
    my $href=$link->attr('href');
    next unless $href =~ s{^\Q$basepath\E}{};
    $href =~ s/#.*$//;  # Remove any #-suffix
    next unless $href;  # Skip if it's empty
    next if exists $unique_links{$href};  # Skip if we've already seen it
    $unique_links{$href} = 1;
    my $page_url = $url . $href;
    my $page_response = $ua->get($page_url);
    
    if ($page_response->is_success) {
        append_md(\$markdown_content, $page_response); # Append to markdown content

    }
}

## Save to a file
#open my $fh, '>', 'output.md';
#print $fh $markdown_content;
#close $fh;
print $markdown_content; # Print to STDOUT

sub append_md {
    my ($markdown_content, $page_response) = @_;
    my $page_tree = HTML::TreeBuilder->new_from_content($page_response->decoded_content);
    my $content = $page_tree->look_down(_tag => 'td', id => 'content');
    if ($content) {
        my $formatter = HTML::FormatMarkdown->new;
        my $markdown = $formatter->format($content);
        $$markdown_content .= $markdown . "\n\n";
    }
    #my $content = $page_tree->look_down(id => 'content');
    #my $page_title = $page_tree->look_down(_tag => 'h1')->as_text;
    $page_tree->delete;
}