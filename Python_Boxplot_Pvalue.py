# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 21:47:47 2020

@author: nitish

Python tool statannot for adding P-value, https://github.com/webermarcolivier/statannot
Integrated statistical tests (binding to scipy.stats methods):
Mann-Whitney
t-test (independent and paired)
Welch's t-test
Levene test
Wilcoxon test
Kruskal-Wallis test
"""
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation

sns.set(style="whitegrid")
df = sns.load_dataset("tips")

x = "day"
y = "total_bill"
hue = "smoker"
ax = sns.boxplot(data=df, x=x, y=y, hue=hue)
ax = sns.stripplot(x=x, y=y, data=df, hue=hue, jitter=0.25, dodge=True)

add_stat_annotation(ax, data=df, x=x, y=y, hue=hue,
                    box_pairs=[(("Thur", "No"), ("Fri", "No")),
                                 (("Sat", "Yes"), ("Sat", "No")),
                                 (("Sun", "No"), ("Thur", "Yes"))
                                ],
                    test='t-test_ind', text_format='full', loc='inside', verbose=2)
"""	
# Use Mann-Whitney test and use star for significant P-value			
add_stat_annotation(ax, data=df, x=x, y=y, hue=hue,
                    box_pairs=[(("Thur", "No"), ("Fri", "No")),
                                 (("Sat", "Yes"), ("Sat", "No")),
                                 (("Sun", "No"), ("Thur", "Yes"))
                                ],
                    test='Mann-Whitney', text_format='star', loc='inside', verbose=2)
"""				
handles, labels = ax.get_legend_handles_labels()

plt.legend(handles[0:2], labels[0:2],loc='upper left', bbox_to_anchor=(1.03, 1))