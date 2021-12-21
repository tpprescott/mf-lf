function make_MF_vertex(node::DecisionTree.Node)
    return MFNode(node.featid, node.featval, make_MF_vertex(node.left), make_MF_vertex(node.right))
end
function make_MF_vertex(leaf::DecisionTree.Leaf)
    return MFLeaf()
end

function get_MF_leaf(μ::MFMean, p::MFIteration)
    return get_MF_leaf(μ, p.θ, p.y_lo)
end
function get_MF_leaf(μ::MFMean, θ, y_lo)
    return get_MF_leaf(μ, vcat(θ, y_lo))
end
function get_MF_leaf(μ::MFMean, features)
    return get_MF_leaf(μ.root, features)
end
function get_MF_leaf(node::MFNode, features)
    child = features[node.featid] < node.featval ? node.left : node.right
    return get_MF_leaf(child, features)
end
function get_MF_leaf(leaf::MFLeaf, features)
    return leaf
end

Base.length(node::MFNode) = length(node.left) + length(node.right)
Base.length(::MFLeaf) = 1

leaves(node::MFNode) = Iterators.flatten((leaves(node.left), leaves(node.right)))
leaves(leaf::MFLeaf)::Tuple{MFLeaf} = (leaf, )
leaves(μ::MFMean) = leaves(μ.root)

function (μ::MFMean)(p)
    l = get_MF_leaf(μ, p)
    return exp(l.logμ)
end
(μ::MFMean)(θ, y_lo) = μ(vcat(θ, y_lo))

function DecisionTree.print_tree(io::IO, leaf::MFLeaf, _k::Int64, depth=-1, indent=0)
    print(io, "\\State \\Return \$\\nu_{$(_k)} = $(round(exp(leaf.logμ), digits=3))\$.\n")
    return _k+1
end
function DecisionTree.print_tree(io::IO, node::MFNode, _k::Int64, depth=-1, indent=0)
    if depth==indent
        print(io, "\n")
        return
    end

    θ_flag = (node.featid <= 3)
    _n = θ_flag ? node.featid : (node.featid - 3)

    if θ_flag
        pars = ["k_{1}", "k_{-1}", "k_2"]
        if_str = "\\If{\$ $(pars[node.featid]) \\leq $(round(node.featval, digits=3)) \$\n}"
    else
        if_str = "\\If{\$ y_{$(_n)} \\leq $(round(node.featval, digits=3)) \$}\n"
    end
    print(io, if_str)
    _k = print_tree(io, node.left, _k, depth, indent+1)
    print(io, "\\Else\n")
    _k = print_tree(io, node.right, _k, depth, indent+1)
    print(io, "\\EndIf\n")
    return _k
end
function DecisionTree.print_tree(io::IO, μ::MFMean, args...)
    print(io, "\\begin{algorithmic}\n")
    print(io, "\\Require \$\\theta=(k_1, k_{-1}, k_2)\$; \$y_{\\lo} = (y_1, y_2, \\dots, y_{10})\$.\n")
    print_tree(io, μ.root, 1, args...)
    print(io, "\\end{algorithmic}\n")
    return true
end
