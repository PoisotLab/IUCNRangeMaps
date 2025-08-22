using DataFrames
import CSV

files = readdir("outputs", join=true)

df = DataFrame.(CSV.File.(files))

full_df = vcat(df...)
CSV.write("bii.csv", full_df)
