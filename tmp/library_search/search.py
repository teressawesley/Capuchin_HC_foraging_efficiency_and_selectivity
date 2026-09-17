from pathlib import Path
import pypdf,json,re
root=Path('C:/Users/TessW/OneDrive/Documents/A - UniKonstanz/Thesis/Thesis Writing/Zotero contents')
out=Path('tmp/library_search')
records=[]
for f in root.rglob('*.pdf'):
 try:
  pages=[p.extract_text() or '' for p in pypdf.PdfReader(f).pages]
  (out/(f.parent.name+'.txt')).write_text('\n\n'.join('PDF PAGE '+str(i+1)+'\n'+t for i,t in enumerate(pages)),encoding='utf-8')
  records.append({'id':f.parent.name,'file':str(f),'pages':len(pages)})
  if re.search(r'stone|hammer|nut.crack',f.name,re.I):
   print('\nFILE',f.parent.name,f.name,flush=True)
   for i,t in enumerate(pages):
    for m in re.finditer(r'success|\d\s*%',t,re.I):
     s=t[max(0,m.start()-140):m.end()+240].replace('\n',' ')
     if 'success' in s.lower(): print('PAGE',i+1,s,flush=True)
 except Exception as e: print('ERROR',f,str(e),flush=True)
(out/'index.json').write_text(json.dumps(records),encoding='utf-8')
