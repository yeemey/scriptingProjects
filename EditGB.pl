#Converts sequence.gb to FASTA
#Uses country, accession, and year info from original gb file as taxa id in FASTA output
#KnownBug1: Country info can be inaccurate when name is subset of longer name.
#For example, Nigeria may be misidentified as Niger.
#KnownBug2: If country uses more than one line in original gb file -> infinite loop(?)
#Current workaround is to quit the script, modify gb file, and rerun.

#!/usr/bin/perl
use strict;
use warnings;

my %CTRY2code=(           ## http://www.iso.org/iso/country_codes/iso_3166_code_lists.htm
 'AFGHANISTAN' => 'AF',
 'ALAND ISLANDS' => 'AX',
 'ALBANIA' => 'AL',
 'ALGERIA' => 'DZ',
 'AMERICAN SAMOA' => 'AS',
 'ANDORRA' => 'AD',
 'ANGOLA' => 'AO',
 'ANGUILLA' => 'AI',
 'ANTARCTICA' => 'AQ',
 'ANTIGUA AND BARBUDA' => 'AG',
 'ARGENTINA' => 'AR',
 'ARMENIA' => 'AM',
 'ARUBA' => 'AW',
 'AUSTRALIA' => 'AU',
 'AUSTRIA' => 'AT',
 'AZERBAIJAN' => 'AZ',
 'BAHAMAS' => 'BS',
 'BAHRAIN' => 'BH',
 'BANGLADESH' => 'BD',
 'BARBADOS' => 'BB',
 'BELARUS' => 'BY',
 'BELGIUM' => 'BE',
 'BELIZE' => 'BZ',
 'BENIN' => 'BJ',
 'BERMUDA' => 'BM',
 'BHUTAN' => 'BT',
 'BOLIVIA' => 'BO',
 'BONAIRE, SINT EUSTATIUS AND SABA' => 'BQ',
 'BOSNIA AND HERZEGOVINA' => 'BA',
 'BOTSWANA' => 'BW',
 'BOUVET ISLAND' => 'BV',
 'BRAZIL' => 'BR',
 'BRUNEI DARUSSALAM' => 'BN',
 'BULGARIA' => 'BG',
 'BURKINA FASO' => 'BF',
 'BURUNDI' => 'BI',
 'CAMBODIA' => 'KH',
 'CAMEROON' => 'CM',
 'CANADA' => 'CA',
 'CAPE VERDE' => 'CV',
 'CAYMAN ISLANDS' => 'KY',
 'CENTRAL AFRICAN REPUBLIC' => 'CF',
 'CHAD' => 'TD',
 'CHILE' => 'CL',
 'CHINA' => 'CN',
 'CHRISTMAS ISLAND' => 'CX',
 'COCOS ISLANDS' => 'CC',
 'KEELING ISLANDS' => 'CC',
 'COLOMBIA' => 'CO',
 'COMOROS' => 'KM',
 'CONGO' => 'CG',
 'THE DEMOCRATIC REPUBLIC OF THE CONGO' => 'CD',
 'COOK ISLANDS' => 'CK',
 'COSTA RICA' => 'CR',
 'COTE D\'IVOIRE' => 'CI',
 'CROATIA' => 'HR',
 'CUBA' => 'CU',
 'CURACAO' => 'CW',
 'CYPRUS' => 'CY',
 'CZECH REPUBLIC' => 'CZ',
 'DENMARK' => 'DK',
 'DJIBOUTI' => 'DJ',
 'DOMINICA' => 'DM',
 'DOMINICAN REPUBLIC' => 'DO',
 'ECUADOR' => 'EC',
 'EGYPT' => 'EG',
 'EL SALVADOR' => 'SV',
 'EQUATORIAL GUINEA' => 'GQ',
 'ERITREA' => 'ER',
 'ESTONIA' => 'EE',
 'ETHIOPIA' => 'ET',
 'FALKLAND ISLANDS' => 'FK',
 'MALVINAS' => 'FK',
 'FAROE ISLANDS' => 'FO',
 'FIJI' => 'FJ',
 'FINLAND' => 'FI',
 'FRANCE' => 'FR',
 'FRENCH GUIANA' => 'GF',
 'FRENCH POLYNESIA' => 'PF',
 'FRENCH SOUTHERN TERRITORIES' => 'TF',
 'GABON' => 'GA',
 'GAMBIA' => 'GM',
 'GEORGIA' => 'GE',
 'GERMANY' => 'DE',
 'GHANA' => 'GH',
 'GIBRALTAR' => 'GI',
 'GREECE' => 'GR',
 'GREENLAND' => 'GL',
 'GRENADA' => 'GD',
 'GUADELOUPE' => 'GP',
 'GUAM' => 'GU',
 'GUATEMALA' => 'GT',
 'GUERNSEY' => 'GG',
 'GUINEA' => 'GN',
 'GUINEA-BISSAU' => 'GW',
 'GUYANA' => 'GY',
 'HAITI' => 'HT',
 'HEARD ISLAND AND MCDONALD ISLANDS' => 'HM',
 'HOLY SEE' => 'VA',
 'VATICAN CITY STATE' => 'VA',
 'HONDURAS' => 'HN',
 'HONG KONG' => 'HK',
 'HUNGARY' => 'HU',
 'ICELAND' => 'IS',
 'INDIA' => 'IN',
 'INDONESIA' => 'ID',
 'ISLAMIC REPUBLIC OF IRAN' => 'IR',
 'IRAQ' => 'IQ',
 'IRELAND' => 'IE',
 'ISLE OF MAN' => 'IM',
 'ISRAEL' => 'IL',
 'ITALY' => 'IT',
 'JAMAICA' => 'JM',
 'JAPAN' => 'JP',
 'JERSEY' => 'JE',
 'JORDAN' => 'JO',
 'KAZAKHSTAN' => 'KZ',
 'KENYA' => 'KE',
 'KIRIBATI' => 'KI',
 'DEMOCRATIC PEOPLE\'S REPUBLIC OF KOREA' => 'KP',
 'REPUBLIC OF KOREA' => 'KR',
 'SOUTH KOREA' => 'KR',
 'KUWAIT' => 'KW',
 'KYRGYZSTAN' => 'KG',
 'LAO PEOPLE\'S DEMOCRATIC REPUBLIC' => 'LA',
 'LAOS' => 'LA',
 'LATVIA' => 'LV',
 'LEBANON' => 'LB',
 'LESOTHO' => 'LS',
 'LIBERIA' => 'LR',
 'LIBYA' => 'LY',
 'LIECHTENSTEIN' => 'LI',
 'LITHUANIA' => 'LT',
 'LUXEMBOURG' => 'LU',
 'MACAO' => 'MO',
 'THE FORMER YUGOSLAV REPUBLIC OF MACEDONIA' => 'MK',
 'MADAGASCAR' => 'MG',
 'MALAWI' => 'MW',
 'MALAYSIA' => 'MY',
 'MALDIVES' => 'MV',
 'MALI' => 'ML',
 'MALTA' => 'MT',
 'MARSHALL ISLANDS' => 'MH',
 'MARTINIQUE' => 'MQ',
 'MAURITANIA' => 'MR',
 'MAURITIUS' => 'MU',
 'MAYOTTE' => 'YT',
 'MEXICO' => 'MX',
 'FEDERATED STATES OF MICRONESIA' => 'FM',
 'REPUBLIC OF MOLDOVA' => 'MD',
 'MONACO' => 'MC',
 'MONGOLIA' => 'MN',
 'MONTENEGRO' => 'ME',
 'MONTSERRAT' => 'MS',
 'MOROCCO' => 'MA',
 'MOZAMBIQUE' => 'MZ',
 'MYANMAR' => 'MM',
 'NAMIBIA' => 'NA',
 'NAURU' => 'NR',
 'NEPAL' => 'NP',
 'NETHERLANDS' => 'NL',
 'NEW CALEDONIA' => 'NC',
 'NEW ZEALAND' => 'NZ',
 'NICARAGUA' => 'NI',
 'NIGER' => 'NE',
 'NIGERIA' => 'NG',
 'NIUE' => 'NU',
 'NORFOLK ISLAND' => 'NF',
 'NORTHERN MARIANA ISLANDS' => 'MP',
 'NORWAY' => 'NO',
 'OMAN' => 'OM',
 'PAKISTAN' => 'PK',
 'PALAU' => 'PW',
 'PALESTINIAN TERRITORY, OCCUPIED' => 'PS',
 'PANAMA' => 'PA',
 'PAPUA NEW GUINEA' => 'PG',
 'PARAGUAY' => 'PY',
 'PERU' => 'PE',
 'PHILIPPINES' => 'PH',
 'PITCAIRN' => 'PN',
 'POLAND' => 'PL',
 'PORTUGAL' => 'PT',
 'PUERTO RICO' => 'PR',
 'QATAR' => 'QA',
 'REUNION' => 'RE',
 'ROMANIA' => 'RO',
 'RUSSIAN FEDERATION' => 'RU',
 'RWANDA' => 'RW',
 'SAINT BARTHELEMY' => 'BL',
 'ST. BARTHELEMY' => 'BL',
 'SAINT HELENA, ASCENSION AND TRISTAN DA CUNHA' => 'SH',
 'ST. HELENA, ASCENSION AND TRISTAN DA CUNHA' => 'SH',
 'SAINT KITTS AND NEVIS' => 'KN',
 'ST. KITTS AND NEVIS' => 'KN',
 'SAINT LUCIA' => 'LC',
 'ST. LUCIA' => 'LC',
 'SAINT MARTIN' => 'MF',
 'ST. MARTIN' => 'MF',
 'SAINT PIERRE AND MIQUELON' => 'PM',
 'ST. PIERRE AND MIQUELON' => 'PM',
 'SAINT VINCENT AND THE GRENADINES' => 'VC',
 'ST. VINCENT AND THE GRENADINES' => 'VC',
 'SAMOA' => 'WS',
 'SAN MARINO' => 'SM',
 'SAO TOME AND PRINCIPE' => 'ST',
 'SAUDI ARABIA' => 'SA',
 'SENEGAL' => 'SN',
 'SERBIA' => 'RS',
 'SEYCHELLES' => 'SC',
 'SIERRA LEONE' => 'SL',
 'SINGAPORE' => 'SG',
 'SINT MAARTEN' => 'SX',
 'SLOVAKIA' => 'SK',
 'SLOVENIA' => 'SI',
 'SOLOMON ISLANDS' => 'SB',
 'SOMALIA' => 'SO',
 'SOUTH AFRICA' => 'ZA',
 'SOUTH GEORGIA AND THE SOUTH SANDWICH ISLANDS' => 'GS',
 'SOUTH SUDAN' => 'SS',
 'SPAIN' => 'ES',
 'SRI LANKA' => 'LK',
 'SUDAN' => 'SD',
 'SURINAME' => 'SR',
 'SVALBARD AND JAN MAYEN' => 'SJ',
 'SWAZILAND' => 'SZ',
 'SWEDEN' => 'SE',
 'SWITZERLAND' => 'CH',
 'SYRIAN ARAB REPUBLIC' => 'SY',
 'TAIWAN' => 'TW',
 'TAJIKISTAN' => 'TJ',
 'UNITED REPUBLIC OF TANZANIA' => 'TZ',
 'THAILAND' => 'TH',
 'TIMOR-LESTE' => 'TL',
 'TOGO' => 'TG',
 'TOKELAU' => 'TK',
 'TONGA' => 'TO',
 'TRINIDAD AND TOBAGO' => 'TT',
 'TUNISIA' => 'TN',
 'TURKEY' => 'TR',
 'TURKMENISTAN' => 'TM',
 'TURKS AND CAICOS ISLANDS' => 'TC',
 'TUVALU' => 'TV',
 'UGANDA' => 'UG',
 'UKRAINE' => 'UA',
 'UNITED ARAB EMIRATES' => 'AE',
 'UNITED KINGDOM' => 'GB',
 'UNITED STATES' => 'US',
 'USA' => 'US',
 'US' => 'US',
 'UNITED STATES MINOR OUTLYING ISLANDS' => 'UM',
 'URUGUAY' => 'UY',
 'UZBEKISTAN' => 'UZ',
 'VANUATU' => 'VU',
 'BOLIVARIAN REPUBLIC OF VENEZUELA' => 'VE',
 'VIET NAM' => 'VN',
 'BRITISH VIRGIN ISLANDS' => 'VG',
 'U.S. VIRGIN ISLANDS' => 'VI',
 'WALLIS AND FUTUNA' => 'WF',
 'WESTERN SAHARA' => 'EH',
 'YEMEN' => 'YE',
 'ZAMBIA' => 'ZM',
 'ZIMBABWE' => 'ZW');
my @CTRYkeys = keys %CTRY2code;                ## Array of country names from the hash
my @FASTAFILE=();                              ## Array of final results
my $FASTAseq=">";                              ## FASTA sequence header result
my $line;                                      ## For reading in GenBank file
my $usergb=$ARGV[0];                           ## User input GenBank file name
chomp $usergb;

open GBFILE,$usergb || die "File not found!";
while ($line = <GBFILE>){
 my $GB2fasta="";                                        ## Extract GB accession number
 my $ctry="";                                            ## Extract country
 my $code="";                                            ## 2-letter country code
 if ($line =~m/ACCESSION/g){                             ## Match to ACCESSION number
  for (my $gbpos=pos($line); substr($line,$gbpos,1) ne "\n"; $gbpos++){
   $GB2fasta =substr($line,$gbpos,1);                    ## Getting GB accession number
   $GB2fasta =~s/\s//g;                                  ## Remove whitespace
   $FASTAseq.=$GB2fasta;                                 ## Concatenate to FASTA header
  }
 }
 elsif ($line =~m/\/country\=\"/g){                      ## Match to country
  for (my $ctrypos=pos($line); substr($line,$ctrypos,1) ne "\""; $ctrypos++){
   $ctry.=substr($line,$ctrypos,1);                      ## Extract country name
  }
  print "$ctry\n";
  my $arraysize = scalar @CTRYkeys;
  my $i=0;
  while ($i<$arraysize){                       ## Search for country name in hash
   if ($ctry =~m/$CTRYkeys[$i]/ig || $CTRYkeys[$i] =~m/$ctry/ig){   
    my $key = $CTRYkeys[$i];                   ## When country matched, get 2-letter code
    $code = $CTRY2code{$key};
    $i=$arraysize+5;
   }
   else {
    $i++;
   }
  }
  $FASTAseq.="_$code";                         ## Append code to FASTA header 
 }
 elsif ($line =~m/\/collection\_date\=/g){     ## Match to collection date
  $line =~m/\d{4}/g;                           ## Match to 4-digit year
  my $yrpos=pos($line);
  my $year = substr($line,$yrpos-4,4);
  $FASTAseq.="_$year";                         ## Append year to FASTA header
 }
 elsif ($line =~m/ORIGIN/g){                   ## Match to sequence 
  my @SEQLINES=();                             ## Array of lines of sequence per taxa
  my $sequenceline;
  $line=<GBFILE>;
  while ($line ne "\/\/\n"){                   ## Before the end of the taxa record...
   $sequenceline = $line;                      ## ... each line of sequence...
   push(@SEQLINES,$sequenceline);              ## ... is pushed into an array
   $line=<GBFILE>;
  }
  my $sequence = join("",@SEQLINES);           ## Sequence array into string
  $sequence =~s/[0-9]//g;                      ## Remove line numbers...
  $sequence =~s/ //g;                          ## ...and spaces
  my $properseq = uc($sequence);               ## Change to upper case letters
  $FASTAseq.="\n$properseq\n";                 ## Append sequence to FASTA header
  push(@FASTAFILE,$FASTAseq);                  ## Add completed sequence record to array
  $FASTAseq=">";                               ## Reset FASTA sequence header
 } 
}
close GBFILE;

my $outputfile="EGB$usergb.fasta";             ## Output to FASTA file
open OUTFILE,">$outputfile" || die "Output file can't be opened!";
print OUTFILE @FASTAFILE;
close OUTFILE;

