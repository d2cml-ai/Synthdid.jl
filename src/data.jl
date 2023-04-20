const cali_path = joinpath(dirname(@__FILE__), "..", "data", "california_prop99.csv");
const quota_path = joinpath(dirname(@__FILE__), "..", "data", "quota.csv");

function california_prop99()
  return CSV.read(cali_path, DataFrame);
end

function quota()
  return CSV.read(quota_path, DataFrame);
end
# data("california_prop99")