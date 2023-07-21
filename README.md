# CSHL_GroupProject2

### Distributed coding of choice, action and engagement across the mouse brain, [Steinmetz,2019](https://www.nature.com/articles/s41586-019-1787-x)

## Project 2: 
Perform multiple regression that predicts firing rate of one cell as functions of firing of other cells. Computer fano factor of each cell and determine if there’s a relationship between the cell's fano factor and its contribution for the target cell. Cycle through cells and average together the betas and fano factors to determine the average network contribution score. 

**Important terms:**
Fano factor - ***squared standard deviation (variance) divided by mean value (when stim introduced, fano factor decreases from ~1.1 to 1***

## Outline: 

- (Replicate Steinmetz kernel regression)
  - Reduced rank regression 
  - Kernel regression 

- Multiple regression 
  - **Clean up data**
    - Get rid of outliers
    - Normalization for each unit across time  
- **Testing: select one area in one animal.**
  - ***Neuron from Neurons***: Randomly choose the first cell in the set (Cell1)
  - Use randomly chosen cells from the rest of the set (N = 1:50) and feed into a multiple regression model to predict the activity of Cell 1 
    - How many cells accurately predict cell 1? 
    - What are the fano factors associated with each cell? 
  - ***Behavior from Neurons***: Replace Cell 1 with behavior (some continuous measure, e.g., eye movement): Repeat the above exercise but replacing Cell 1 with behavioral measure. 
    - How many cells accurately predict Behavior? 
    - What are the fano factors associated with each cell?
    - Is there a continuous behavior that better predicts activity?  
  - ***How well do the fano factors correlate with the predictive values?***
- **Scaling up: Apply to multiple brain areas** (if we have time) 
  - Reduce neurons in all areas (or a handful, e.g., N = 1:20) into a single estimate of neural activity (so 1 estimate per area)
  - Use multiple regression to determine if a collection of areas accurately predict behavior (use the best behavior  from testing) 
    - Iterate through this to have different combinations of areas
      
## MXC Suggestions
- In the regression, penalize regressors 
- Formal model comparison (penalize regressions with more regressors) 
- Cross-validation
- Draw out ideas for final figures 
