# Pandas

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Read data](#read-data)
	- [From excel sheets](#from-excel-sheets)
- [Write data](#write-data)
	- [To CSV](#to-csv)
		- [Just one column](#just-one-column)
	- [To Excel sheets](#to-excel-sheets)
- [Data types](#data-types)
- [indexes](#indexes)
	- [Set indexes](#set-indexes)
	- [Set multiple indexes](#set-multiple-indexes)
	- [drop (reset) indexes](#drop-reset-indexes)
- [Columns and rows](#columns-and-rows)
	- [Columns names](#columns-names)
		- [Assign names (all columns)](#assign-names-all-columns)
		- [Renaming](#renaming)
		- [Get column names as lists](#get-column-names-as-lists)
			- [Simple dataframe](#simple-dataframe)
			- [Multi indexed columns](#multi-indexed-columns)
	- [Select columns and rows](#select-columns-and-rows)
		- [With indexes](#with-indexes)
			- [General case](#general-case)
			- [Select one value](#select-one-value)
		- [Columns](#columns)
			- [Select one column](#select-one-column)
				- [as data frame](#as-data-frame)
				- [as series](#as-series)
				- [as list](#as-list)
			- [Select several columns as data frame](#select-several-columns-as-data-frame)
			- [Drop Columns](#drop-columns)
		- [Select rows](#select-rows)
			- [By slice (and specific columns)](#by-slice-and-specific-columns)
			- [By numeric values](#by-numeric-values)
				- [Operators: <, >, ==, !=, <=, >=](#operators-)
				- [Operators: & (and), | (or), ~ (not)](#operators-and-or-not)
			- [By regular expression](#by-regular-expression)
			- [Select with a list](#select-with-a-list)
				- [Is in the list](#is-in-the-list)
				- [Is not in the list](#is-not-in-the-list)
	- [Operations on Columns](#operations-on-columns)
	- [Unique values](#unique-values)
		- [Unique values from a column into a list](#unique-values-from-a-column-into-a-list)
		- [Number of values per category](#number-of-values-per-category)
		- [Count the number of unique values (categories)](#count-the-number-of-unique-values-categories)
- [Strings](#strings)
	- [All uppercase](#all-uppercase)
	- [Replace characters](#replace-characters)
	- [Split and add to new Columns](#split-and-add-to-new-columns)
	- [Fill numericals with zeros](#fill-numericals-with-zeros)
	- [Concatenate integers as Strings](#concatenate-integers-as-strings)

<!-- /TOC -->

## Read data
### From excel sheets
```python
samples = pd.read_excel('~/path/to/SequenceGroups.xlsx',sheet_name='List1ForPaper')
```
## Write data
### To CSV

#### Just one column
* sometimes needs quoting=csv.QUOTE_NONE. Sometimes not!
```python
df['Fasta'].to_csv('/Users/avignal/GenotoulWork/BlastSNP_HAv3/Fasta.csv', index=False, header=False, escapechar=' ', quoting=csv.QUOTE_NONE)
```

### To Excel sheets
```python
path = '/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/ListsSamples/SequenceGroups2.xlsx'
writer = pd.ExcelWriter(path, engine='xlsxwriter')
xlsx = pd.ExcelFile(path)
for sheet in xlsx.sheet_names:
    df = xlsx.parse(sheet_name=sheet, index_col=0)
    df.to_excel(writer, sheet_name=sheet)
eigenvec.to_excel(writer, sheet_name='Vec')
eigenval.to_excel(writer, sheet_name='Val')
writer.save()
writer.close()
```

## Data types
* of a dataframe columns

```python
df.dtypes
```

## indexes
### Set indexes
```python
assembly = assembly.set_index('GenbankName')
```
### Set multiple indexes
```python
blast370index3 = blast370.set_index(keys=['Query','Subject','sPos','qStart']).sort_index(ascending=[True,True,True,True])
```
### drop (reset) indexes
```python
df.reset_index()
```

## Columns and rows
### Columns names
#### Assign names (all columns)
```python
eigenvec.columns = ['name','name2','PC1','PC2','PC3']
```
#### Renaming
* Uses a dictionary: {odlName1:newName1,odlName2:newName2...odlNameN:newNameN}

```python
Pig150F11.rename(columns={'Hits':'Pig150F11_hits'})
```

#### Get column names as lists
##### Simple dataframe

```python
list(selStats2.columns)
```

##### Multi indexed columns

```python
list(selStats2.columns.levels[1])
```

### Select columns and rows

#### With indexes
##### General case
```python
df.iloc[rows, columns]
df.loc[rows, columns]
```

##### Select one value
* To assign to a variable
```Python
eigenvec = eigenvec.set_index('name')
var = eigenvec.at['YC9','PC2']
```
* To assign a new values
```Python
eigenvec = eigenvec.set_index('name')
eigenvec.at['YC9','PC2'] = 'new value'
```

#### Columns
##### Select one column
###### as data frame
```python
samples.loc[:,['IndSeqName']]
```

###### as series
```python
samples['IndSeqName']
```

###### as list
```python
df10.loc[:,['fasta']].tolist()
```

##### Select several columns as data frame

```python
samples.loc[:,['IndSeqName','Ref','Type','Loc','Country','Type.1','Type.2','Region']]
```

* With column indexes
```python
college.iloc[:, [4,6]].head()
```
* With column Names
```python
college.loc[:, [‘WOMENONLY’, ‘SATVRMID’]].head()
```

##### Drop Columns
```python
eigenvec = eigenvec.drop(columns=['name2','col2'])
```

#### Select rows
##### By slice (and specific columns)
* By column names
```python
data.loc[1:2,['One','Four']]
```
* By column index
```python
data.iloc[1:3,[4,8]]
```

##### By numeric values
###### Operators: <, >, ==, !=, <=, >=
```python
df2 = df1[(df1.Period_Size > 10)]
```

###### Operators: & (and), | (or), ~ (not)
*Take precedence over the comparison operators. Need parentheses.

* One line
```python
df2 = df1[(df1.Copy_Number > 50) & (df1.Period_Size > 80) & (df1.Period_Size < 110)]
```

* More lisible
```python
Criteria1 = df1.Copy_Number > 50
Criteria1 = df1.Period_Size > 80
Criteria2 = df1.Period_Size < 110
Criteria = Criteria1 & Criteria2 & Criteria"
df2 = df1[Criteria]
```

##### By regular expression
* Select lines of the eigenvec table, with values in the names column starting with different patterns
```python
eigenvec.loc[eigenvec.name.str.contains('^FL|^NM|^XC|^YC',regex=True),:]
```
* Assign a colour these entries
```python
eigenvec.loc[eigenvec.name.str.contains('^FL|^NM|^XC|^YC',regex=True),'Color'] = '#ecec13ff'
```

##### Select with a list
###### Is in the list
```Python
telomerePentas = ["GGTTA","GTTAG","TTAGG","TAGGT","AGGTT","TAACC","AACCT","ACCTA","CCTAA","CTAAC"]
df2 = df1[df1['Repseq'].isin(telomerePentas)]
```

###### Is not in the list
```python
telomerePentas = ["GGTTA","GTTAG","TTAGG","TAGGT","AGGTT","TAACC","AACCT","ACCTA","CCTAA","CTAAC"]
df2 = df1[~df1['Repseq'].isin(telomerePentas)]
```

### Operations on Columns
* Will create a third column if non-existent

```python
HSAchrs['Pig150F11_perMb'] = HSAchrs['Pig150F11_hits'] / HSAchrs['Mb']
```

### Unique values
#### Unique values from a column into a list
```python
chromosomes = df1['Chr'].unique().tolist()
```
#### Number of values per category
```python
dataSamples['CategPCA'].value_counts()
```

#### Count the number of unique values (categories)
* Number of categories
```python
dataSamples['CategPCA'].nunique()
```

## Strings
### All uppercase
```python
nz1['Ruche'] = nz1.hive_number_bs.str.upper()
```

### Replace characters
```python
nz1['Ruche'] = nz1.Ruche.str.replace( "_", "-", n=-1, case=None, flags=0, regex=True)
```

### Split and add to new Columns
* n = the number of splits to perform
* The default expand=True will return a series of lists
* expand=True will return a dataframe
```python
sequences[['Year','Number']] = sequences['Ruche'].str.split("-", n=1, expand=True)
```
* If only one part of the string is needed
```python
samples['tag'] = samples['IndSeqName'].str.split("-", n=2, expand=True)[1]
```

### Fill numericals with zeros
```python
sequences['Number'] = sequences['Number'].str.zfill(4)
```

### Concatenate integers as Strings

```python3
plinkData['forJoin'] = plinkData['chrPlink'].astype(str) + ":" + plinkData['pos'].astype(str)
```



end of file
