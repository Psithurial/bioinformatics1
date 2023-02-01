#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 11:41:18 2022

@author: junevolkmann
"""

#interface
import random

#file read
import os.path
fPath = "/home/junevolkmann/Downloads/"
address = fPath + "dataset_163_4.txt"
exists = os.path.exists(address)
if exists:
  file = open(address, 'r', encoding='utf-8-sig')
  Line = []
  Line = file.readlines()

#README
#text = string
#pattern = string < text
#k = len(pattern)
#n = len(text)
#li = list
#l = 
#d = hamming distance
#t = len(text_list)

#neucleotides
nT = ['A','T','C','G']

#common variables:
#returns all possible beginning
#k-mer indices in a genome length n)
def mer_ind(n, k):
  return range(0, n - k + 1)

#returns all possible mers of length k ias strings in a list
def all_mers(k):
  mer = ""
  for i in range(0, k):
    mer = mer + "A"
  return neighbors(mer, k)

#returns a window of k characters in string t from  index i
def window(text, i, k):
  return text[i:i + k]
#end common variables

#data formatting

#returns a list of strings from a text file
def remove_newlines(fname):
  flist = open(fname).readlines()
  return [s.rstrip('\n') for s in flist]

def list_to_string(li):
  st = ""
  for i in li:
    st = st + str(i) + " "
  return st

def string_to_list(st):
  li = st.split(" ")
  return li

def list_from_file_lines(start, end):
  li = []
  for i in range(start, end):
    li.append(Line[i].rstrip('\n'))
  return li

def read_profile_matrix_from_file(line_num):
  matrix = {"A": [], "C": [], "T": [], "G": []}
  matrix["A"] = (Line[line_num].rstrip('\n').split(" "))
  matrix["C"] = (Line[line_num + 1].rstrip('\n').split(" "))
  matrix["G"] = (Line[line_num + 2].rstrip('\n').split(" "))
  matrix["T"] = (Line[line_num + 3].rstrip('\n').split(" "))
  return matrix
  
#end data formatting
#week 1
def frequency_map(text, k):
  te = text.strip()
  n = len(text)
  k = int(k)
  r = mer_ind(n, k)
  f_map = {}
  for i in r:
      p = window(te, i, k)
      f_map[p] = 0
      for j in r:
        w = window(te, j, k)
        if w == p:
          f_map[p] = f_map[p] + 1
  return f_map


def frequent_words(text, k):
  te = text.strip()
  k = int(k)
  words = []
  f = frequency_map(te, k)
  m = max(f.values())
  for i in f:
      if f[i] == m:
        words.append(i)
  return words
  
def reverse(pattern):
  rev = ""
  p = pattern.strip()
  for i in p:
    rev = i + rev
  return rev

def complement(pattern):
  com = ""
  p = pattern.strip()
  for i in p:
    if i == "A":
      com = com + "T"
    elif i == "T":
      com = com + "A"
    elif i == "C":
      com = com + "G"
    elif i == "G":
      com = com + "C"
  return com
  
def reverse_complement(pattern):
    p = pattern.strip()
    p = reverse(p)
    p = complement(p)
    return p
  
def pattern_matching(pattern, text):
  p = pattern.strip()
  te = text.strip()
  k = len(p)
  n = len(te)
  r = mer_ind(n, k)
  pos = []
  for i in r:
      w = window(te, i, k)
      if w == p:
          pos.append(i)
  return pos
  
def pattern_count(pattern, text):
  p = pattern.strip()
  te = text.strip()
  k = len(p)
  n = len(te)
  r = mer_ind(n, k)
  count = 0
  for i in r:
    w = window(te, i, k)
    if w == p:
      count = count + 1
  return count

def max_map(freq_map):
  m = max(freq_map.values())
  return m

def better_frequent_words(text, k):
  te = text.strip()
  k = int(k)
  f_pats = []
  f_map = frequency_map(te, k)
  m = max_map(f_map)
  for pat in f_map:
    if f_map[pat] == m:
      f_pats.append(pat)
  return f_pats
  
def find_clumps(text, k, l, d):
  te = text.strip()
  n = len(te)
  k = int(k)
  l = int(l)
  d = int(d)
  pats = []
  r = mer_ind(n, l)
  for i in r:
    w = window(te, i, l)
    f_map = frequency_map(w, k)
    for s in f_map:
      if f_map[s] >= d:
        pats.append(s)
  pats = [*set(pats)]
  return pats


#week 2

def count_mers(text, k):
  te = text.strip()
  n = len(te)
  k = int(k)
  r = mer_ind(n, k)
  mers = []
  counts = {}
  for i in r:
    w = window(te, i, k)
    if w in mers:
      counts[w] = counts[w] + 1
    else:
      mers.append(w)
      counts[w] = 0
  mode = (max(counts, key=counts.get))
  return counts, mode

def hamming(text1, text2):
  t1 = text1.strip()
  t2 = text2.strip()
  n = len(t1)
  k = len(t2)
  d = 0
  if n != k:
    print("Error: hamming params len !=")
    print(str(n) + " " + t1)
    print(str(k) + " " + t2)
  else:
    for i in range(0, n):
      if t1[i] != t2[i]:
        d = d + 1
    return d

def skew(text, char1, char2):
  t = text.strip()
  c1 = char1.strip()
  c2 = char2.strip()
  n = len(t)
  sN = 0
  sI = [0]
  sMin = 0
  sMins = []
  sMax = 0
  sMaxs = []
  for i in range(0, n):
    if t[i] == c1:
      sN = sN + 1
    elif t[i] == c2:
      sN = sN - 1
    sI.append(sN)
  for i in range(0, len(sI)):
    if sI[i] <= sMin:
      sMin = sI[i]
    elif sI[i] >= sMax:
      sMax = sI[i]
  for i in range(0, len(sI)):
    if sI[i] == sMin:
      sMins.append(i)
    elif sI[i] == sMax:
      sMaxs.append(i)
  print("Skew i: " + list_to_string(sI))
  print("Skew Min: " + str(sMin))
  print("Skew Min Index: " + str(sMins))
  print("Skew Max: " + str(sMax))
  print("Skew Max Index:: " + str(sMaxs))

def approximate_match_locations(pattern, text, d):
  p = pattern.strip()
  te = text.strip()
  k = len(p)
  n = len(te)
  d = int(d)
  r = mer_ind(n, k)
  count = 0
  ind = []
  for i in r:
    w = window(te, i, k)
    if hamming(w, p) <= d:
      count = count + 1
      ind.append(i)
  return ind

def approximate_match_count(pattern, text, d):
  p = pattern.strip()
  te = text.strip()
  k = len(p)
  n = len(te)
  d = int(d)
  r = mer_ind(n, k)
  count = 0
  ind = []
  for i in r:
    w = window(te, i, k)
    if hamming(w, p) <= d:
      count = count + 1
      ind.append(i)
  return count

def motif_counts(text, k, d):
  te = text.strip()
  n = len(te)
  k = int(k)
  d = int(d)
  r = mer_ind(n, k)
  motifs = []
  m_counts = {}
  motifs.append(text[0:k])
  m_counts[motifs[0]] = 1
  for i in r:
    w = window(te, i, k)
    v = range(0, len(motifs))
    matches = 0
    for j in v:
      m = motifs[j]
      if hamming(w, m) <= d:
          matches = matches + 1
          m_counts[m] = m_counts[m] + 1
    if matches == 0:
        motifs.append(w)
        m_counts[w] = 1
  return m_counts

def suffix(text):
    s = text[1:]
    return(s)

def neighbors(text, d):
  te = text.strip()
  n = len(te)
  d = int(d)
  if d == 0:
    return te
  if n == 1:
    return nT
  s = suffix(te)
  p = text[0]
  h = []
  suf_h = neighbors(s, d)
  for i in suf_h:
    if hamming(s, i) < d:
      for j in range(0, len(nT)):
        h.append(nT[j] + i)
    else:
      h.append(p + i)
  return h

def d_neighbors(text, d):
  te = text.strip()
  n = len(text)
  d = int(d)
  if d == 0:
    return te
  if n == 1:
    return nT
  s = suffix(te)
  p = te[0]
  h = []
  suf_h = d_neighbors(s, d)
  for i in suf_h:
    if hamming(s, i) == d - 1:
      for j in nT:
        if j != p:
          h.append(j + i)
    elif hamming(s, i) == d:
        h.append(p + i)
  return h

def frequent_words_with_mismatches(text, k, d):
  te = text.strip()
  n = len(te)
  k = int(k)
  d = int(d)
  r = mer_ind(n, k)
  pats = []
  f_map = {}
  for i in r:
    w = window(te, i, k)
    h = neighbors(w, d)
    u = len(h)
    for j in range(0, u):
      g = h[j]
      if g not in f_map:
        f_map[g] = 1
      else:
        f_map[g] = f_map[g] + 1
  m = max_map(f_map)
  for  p in f_map:
    if f_map[p] == m:
            pats.append(p)
  return pats

def frequent_words_w_mismatches_reverse_complements(text, k, d):
  te = text.strip()
  n = len(te)
  k = int(k)
  d = int(d)
  r = mer_ind(n, k)
  pats = []
  f_map = {}
  for i in r:
    w = window(te, i, k)
    wRC = reverse_complement(w)
    h = neighbors(w, d)
    hRC = neighbors(wRC, d)
    u = len(h)
    for j in range(0, u):
      g = h[j]
      b = hRC[j]
      if g not in f_map:
        f_map[g] = 1
      else:
        f_map[g] = f_map[g] + 1
      if b not in f_map:
        f_map[b] = 1
      else:
        f_map[b] = f_map[b] + 1
  m = max_map(f_map)
  for  p in f_map:
    if f_map[p] == m:
            pats.append(p)
  return pats

# Week 3
  
def motif_enumeration(text_list, k, d):
  li = text_list
  k = int(k)
  d = int(d)
  t1 = text_list[1]
  n1 = len(t1)
  pats = []
  for i in mer_ind(n1, k):
    w = window(t1, i, k)
    h = neighbors(w, d)
    for j in h:
      b = neighbors(j, d)
      count_b = 0
      for s in li[0:]:
        for v in b:
          if v in s:
            count_b = count_b + 1
            break
            #if all strings are positive for one or more neighbor of j, j goes in pats
      if count_b == len(li):
        pats.append(j)
  pats = [*set(pats)]
  return pats

#week 3 base
def motifs(text_list, k):
  motifs = []
  k = int(k)
  for l in range(0, len(text_list)):
    motifs.append(text_list[l][0:k])
  return motifs


def count(motifs):
  t = len(motifs)
  k = len(motifs[0])
  count_matrix = {"A" : [], "C" : [], "G" : [], "T" : []}
  for l in count_matrix:
    for i in range(0, k):
      rc = 1
      for j in range(0, t):
        a = motifs[j][i]
        if a == l:
          rc = rc + 1
      count_matrix[l].append(rc)
  return count_matrix

def profile(motifs):
  c = count(motifs)
  n = len(motifs)
  profile = {"A" : [], "C" : [], "G" : [], "T" : []}
  for l in c:
    for i in range(0, len(c[l])):
      profile[l].append((c[l][i]/n/1))
  return profile

def consensus(motifs):
  p = profile(motifs)
  consensus = ""
  nl = "A"
  n = len(p[nl])
  for i in range(0, n):
    for l in p:
      if p[l][i] > p[nl][i]:
        nl = l
    consensus = consensus + nl
  return consensus

def score(motifs):
  score = 0
  t = len(motifs)
  cons = consensus(motifs)
  for i in range(0, t):
    d = hamming(cons, motifs[i])
    score = score + d
  return score
  


def probability_per_profile(pattern, profile):
  prob = 1
  for i in range(0, len(pattern)):
    prob = prob * profile[pattern[i]][i]
  return prob
#end week 3 base
    
def min_motif_score(text, pattern):
  te = text
  p = pattern
  n = len(te)
  k = len(p)
  min_score = k
  for i in mer_ind(n, k):
    w = window(te, i, k)
    d = hamming(w, p)
    if d < min_score:
      min_score = d
  return min_score

def motif_score(text_list, pattern):
  p = pattern
  total_list_score = 0
  for t in text_list:
    s = min_motif_score(t, p)
    total_list_score = total_list_score + s
  return total_list_score

#def motif_matrix(text_list, k):
  #matrix = [text_list[0, len(text_list)][0, k]]
  #return matrix

def median_string(text_list, k):
  n = len(text_list)
  dist = k*n
  median = ""
  for m in all_mers(k):
    d = motif_score(text_list, m)
    if dist > d:
      dist = d
      median = m
  return median

# should return m only when d(m) matches lowest dist
def all_median_strings(text_list, k):
  n = len(text_list)
  dist = k*n
  medians = []
  for m in all_mers(k):
    d = motif_score(text_list, m)
    if dist > d:
      dist = d
      medians.append(m)
  for med in medians:
    d2 = motif_score(text_list, med)
    if d2 > dist:
      del med
  return medians


def profile_most_probable(text, k, profile_matrix):
  te = text
  n = len(te)
  k = int(k)
  max_p = 0
  pmp = window(te, 0, k)
  for i in mer_ind(n, k):
    w = window(te, i , k)
    p = 1
    for j in range(0, k):
      c = w[j]
      mc = float(profile_matrix[c][j])
      p = p * mc
    if max_p < p:
      max_p = p
      pmp = w
  return pmp

def greedy_motif_search(text_list, k, t):
  best_motifs = motifs(text_list, k)
  n = len(text_list[0])
  st1 = text_list[0]
  sts = text_list[1:]
  for i in mer_ind(n, k):
    w = window(st1, i, k)
    mot = [w]
    for j in sts:
      profile_matrix = profile(mot, k)
      next_motif = profile_most_probable(j, k, profile_matrix)
      mot.append(next_motif)
    if score(mot) < score(best_motifs):
      best_motifs = mot
  return best_motifs
#end week 3
#week 4
  
def motifs_profile(profile, text_list):
  k = len(profile["A"])
  t = len(text_list)
  motifs = []
  for i in range(0, t):
    pr = profile_most_probable(text_list[i], k, profile)
    motifs.append(pr)
  return motifs

def randomized_motif_search(text_list, k, t):
  motifs = []
  n = len(text_list[0])
  for l in range(0, t):
    i = random.randint(0, n - k)
    motifs.append(window(text_list[l], i, k))
  best_motifs = motifs
  while True:
    prof = profile(motifs, k)
    motifs = motifs_profile(prof, text_list)
    if score(motifs) < score(best_motifs):
      best_motifs = motifs[:]
    else:
      return best_motifs

def random_motif_search_1k(text_list, k, t):
  best = randomized_motif_search(text_list, k, t)
  for i in range(0, 1000):
    new = randomized_motif_search(text_list, k, t)
    if score(best) > score(new):
      best = new
  return best

def gibbs_sampler(text_list, k, t, N):
  best_motifs = []
  motifs = []
  n = len(text_list[0])
  r = random.randint(0, n - k)
  for l in range(0, len(text_list)):
    motifs.append(window(text_list[l], r, k))
  best_motifs = motifs[:]
  for j in range(0, N):
    i = random.randint(0, t - 1)
    nm = motifs[:]
    del nm[i]
    prof = profile(nm)
    motifs[i] = profile_most_probable(text_list[i], k, prof)
    if score(motifs) < score(best_motifs):
      best_motifs = motifs[:]
  return best_motifs

def gibbs_motifs_list(text_list, k, t, N, j):
  mli = []
  for i in range(0, j):
    mli.append(gibbs_sampler(text_list, k, t, N))
  return mli


def best_motifs(motifs_list):
  mli = motifs_list
  bs = score(mli[0]) #best score
  bm = [] #best motif
  for i in mli:
    s = score(i)
    if s < bs:
      bs = s
      bm = i
  print(bs)
  return bm
 
text_string = "CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
DNA = string_to_list(Line[1].rstrip('\n'))
DNAs = string_to_list(text_string)
#print(list_to_string(best_motifs(gibbs_motifs_list(DNA, 15, 20, 2000, 20))))

print(all_median_strings(["CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC", "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC", "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"], 7))
pro = .3 * .3 * .1 * .5 * .9
print(pro)