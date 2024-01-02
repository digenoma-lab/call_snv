
open(FILE,$ARGV[0]) or die "cannot open file\n";
my $size=5000000;
#hg38
%chrs=(chr1=>1,chr2=>1,chr3=>1,chr4=>1,chr5=>1,chr6=>1,chr7=>1,chr8=>1,chr9=>1,chr10=>1,chr11=>1,chr12=>1,chr13=>1,chr14=>1,chr15=>1,chr16=>1,chr17=>1,chr18=>1,chr19=>1,chr20=>1,chr21=>1,chr22=>1,chrX=>1,chrY=>1,chrM=>1);

while(my $line=<FILE>){
	chomp $line;
	my ($chr, $len,undef)=split("\t",$line);
	next if(!defined $chrs{$chr});

	if($len > $size){
		my $w=int($len/$size);	
		my $j=1;
		for(my $i=1; $i <=$w ; $i++){
			print join("\t",$chr,$j,$j+$size)."\n";
			$j+=$size;
		}
		print join("\t",$chr,$j,$len)."\n";
		
	}else{
		print join("\t",$chr,1,$len)."\n";
	}
}

