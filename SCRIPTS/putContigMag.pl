#
#   pragmas
#
use strict;
use warnings;
my ($id,$seq,$line,@line,%hContig,$magnitud,$PREFIX);
open (SAM, "$ARGV[0]") || die "No se puede abrir el archivo SAM $ARGV[0]" ;
open (CONTIG, "$ARGV[1]") || die  "No se puede abrir el archivo FASTA $ARGV[1]"; 
open (OUTC, ">$ARGV[2]") || die  "No se puede abrir el archivo $ARGV[2]"; 
$PREFIX=$ARGV[3];

print "Entro putContig.pl\n";

print "Processing SAM....\n";
while ($id=<SAM>){
     if($id =~ /^$PREFIX/){
		@line=split(/\t/,$id);
		my $id2 =$line[2];
		$id2 =~ s/\s/_/g;
		if(defined($hContig{$id2})){
			#print "$id2\n";
			$hContig{$id2}[0]=$hContig{$id2}[0]+1;
	    }else{
			$hContig{$id2}[0]=1;
			#print "Inicializar: $id2\n";
		}
	 }
}
print "End\n";
close(SAM);


print "Processing CONTIG....\n";
while ($id=<CONTIG>){
	chomp($id);
	#print "$id\n";
    if($id=~ /^\>contig/ || $id=~ /^\>NODE/){
	    chomp($id);
		$seq=<CONTIG>;
		$id =~ s/(?<=..)\s/_/g;
		$id=substr $id, 1;
		#print "Id a buscar: $id\n";
		if( (defined($hContig{$id}))){
			print OUTC ">$id magnitude=$hContig{$id}[0]\n$seq";
		}
	 }
}
close(OUTC);
close(OUTC)






