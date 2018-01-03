#!/usr/bin/perl -w
$vf = "./version";
if (-f $vf) {
    chomp($version = `cat ./version`);
    $version += 1;
    `/bin/rm -f ./version`;
} else {
    chomp($version = `date +%F%T`);
    $version =~ s/-/_/g;
    $version =~ s/:/_/g;
}
`echo $version > ./version`;
`git add ./version`;
`git commit -m "version change $version"`;
`git archive --format=tar -o ../boltzmann_$version.tar master`;
`gzip ../boltzmann_$version.tar`;
exit(0);
