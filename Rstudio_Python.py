import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import AdaBoostRegressor
import os
path='C:/Users/nitis/Desktop'
os.chdir(path=path)
# Create the dataset
rng = np.random.RandomState(1)
X = np.linspace(0, 6, 100)[:, np.newaxis]
y = np.sin(X).ravel() + np.sin(6 * X).ravel() + rng.normal(0, 0.1, X.shape[0])

# Fit regression model
regr_1 = DecisionTreeRegressor(max_depth=4)
regr_2 = AdaBoostRegressor(DecisionTreeRegressor(max_depth=4),
                          n_estimators=300, random_state=rng)

regr_1.fit(X, y)
regr_2.fit(X, y)

# Predict
y_1 = regr_1.predict(X)
y_2 = regr_2.predict(X)

# Plot the results
plt.figure()
plt.scatter(X, y, c="k", label="training samples")
plt.plot(X, y_1, c="g", label="n_estimators=1", linewidth=2)
plt.plot(X, y_2, c="r", label="n_estimators=300", linewidth=2)
plt.xlabel("data")
plt.ylabel("target")
plt.title("Boosted Decision Tree Regression")
plt.legend()
plt.show()

###################################
from sksurv.datasets import load_veterans_lung_cancer
data_x, data_y = load_veterans_lung_cancer()
#data_y
pd.DataFrame.from_records(data_y[[11, 5, 32, 13, 23]], index=range(1, 6))


from sksurv.nonparametric import kaplan_meier_estimator
time, survival_prob = kaplan_meier_estimator(data_y["Status"], data_y["Survival_in_days"])
plt.step(time, survival_prob, where="post")
plt.ylabel("est. probability of survival $\hat{S}(t)$")
plt.xlabel("time $t$")
plt.show()

data_x["Treatment"].value_counts()
for treatment_type in ("standard", "test"):
    mask_treat = data_x["Treatment"] == treatment_type
    time_treatment, survival_prob_treatment = kaplan_meier_estimator(
        data_y["Status"][mask_treat],
        data_y["Survival_in_days"][mask_treat])

    plt.step(time_treatment, survival_prob_treatment, where="post",
             label="Treatment = %s" % treatment_type)

plt.ylabel("est. probability of survival $\hat{S}(t)$")
plt.xlabel("time $t$")
plt.legend(loc="best")
plt.show()

######################################
from lifelines.datasets import load_dd
data = load_dd()
data.head()

from lifelines import KaplanMeierFitter
kmf = KaplanMeierFitter()

T = data["duration"]
E = data["observed"]
kmf.fit(T, event_observed=E)

from matplotlib import pyplot as plt
kmf.survival_function_.plot()
plt.title('Survival function of political regimes')
plt.show()

from lifelines.utils import median_survival_times
median_ci = median_survival_times(kmf.confidence_interval_)
ax = plt.subplot(111)
dem = (data["democracy"] == "Democracy")
kmf.fit(T[dem], event_observed=E[dem], label="Democratic Regimes")
kmf.plot(ax=ax)
kmf.fit(T[~dem], event_observed=E[~dem], label="Non-democratic Regimes")
kmf.plot(ax=ax)
plt.title("Lifespans of different global regimes")
plt.show()


ax = plt.subplot(111)
t = np.linspace(0, 50, 51)
kmf.fit(T[dem], event_observed=E[dem], timeline=t, label="Democratic Regimes")
ax = kmf.plot(ax=ax)
kmf.fit(T[~dem], event_observed=E[~dem], timeline=t, label="Non-democratic Regimes")
ax = kmf.plot(ax=ax)
plt.title("Lifespans of different global regimes")
plt.show()


from lifelines.statistics import logrank_test
results = logrank_test(T[dem], T[~dem], E[dem], E[~dem], alpha=.99)
results.print_summary()

regime_types = data['regime'].unique()
for i, regime_type in enumerate(regime_types):
    ax = plt.subplot(2, 3, i + 1)
    ix = data['regime'] == regime_type
    kmf.fit(T[ix], E[ix], label=regime_type)
    kmf.plot(ax=ax, legend=False)
    plt.title(regime_type)
    plt.xlim(0, 50)
    if i==0:
        plt.ylabel('Frac. in power after $n$ years')
plt.tight_layout()
plt.show()


kmf = KaplanMeierFitter().fit(T, E, label="all_regimes")
kmf.plot(at_risk_counts=True)
#plt.show()
plt.savefig('kmf.png')## Save images in Python matplotlib


#########################################################
############ Python box-plot with p-value ###############
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
plt.show()



### Save SNS figure in Python #####
import seaborn as sns
sns.set_theme(style="whitegrid")

# Load the example tips dataset
tips = sns.load_dataset("tips")

# Draw a ested violinplot and split the violins for easier comparison
sns_plot = sns.violinplot(data=tips, x="day", y="total_bill", hue="smoker",
               split=True, inner="quart", linewidth=1,
               palette={"Yes": "b", "No": ".85"})
#plt.show()
fig = sns_plot.get_figure()
fig.savefig('SNS_Violin.png')


import pandas as pd
import matplotlib.pyplot as plt
 
# Create a dataframe
value1=np.random.uniform(size=20)
value2=value1+np.random.uniform(size=20)/4
df = pd.DataFrame({'group':list(map(chr, range(65, 85))), 'value1':value1 , 'value2':value2 })
 
# Reorder it following the values of the first value:
ordered_df = df.sort_values(by='value1')
my_range=range(1,len(df.index)+1)
 
# The vertical plot is made using the hline function
# I load the seaborn library only to benefit the nice looking feature
import seaborn as sns
plt.hlines(y=my_range, xmin=ordered_df['value1'], xmax=ordered_df['value2'], color='grey', alpha=0.4)
plt.scatter(ordered_df['value1'], my_range, color='skyblue', alpha=1, label='value1')
plt.scatter(ordered_df['value2'], my_range, color='green', alpha=0.4 , label='value2')
plt.legend()
 
# Add title and axis names
plt.yticks(my_range, ordered_df['group'])
plt.title("Comparison of the value 1 and the value 2", loc='left')
plt.xlabel('Value of the variables')
plt.ylabel('Group')
plt.show()

###################################
############# RPKM ################
from bioinfokit.analys import norm, get_data
### load sugarcane RNA-seq expression dataset (Published in Bedre et al., 2019)
df = get_data('sc_exp').data
df = df.set_index('gene')
nm = norm()
nm.rpkm(df=df, gl='length')
#Get RPKM normalized dataframe
rpkm_df = nm.rpkm_norm
rpkm_df.head()
df.to_csv('bioinfokit.csv')

############## TPM ################
from bioinfokit.analys import norm, get_data
# load sugarcane RNA-seq expression dataset (Published in Bedre et al., 2019)
df = get_data('sc_exp').data
df = df.set_index('gene')
df.head()

nm = norm()
nm.tpm(df=df, gl='length')
# get TPM normalized dataframe
tpm_df = nm.tpm_norm
tpm_df.head(10)
