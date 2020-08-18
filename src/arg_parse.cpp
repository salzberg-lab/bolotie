#include <iostream>
#include <unistd.h>
#include <iomanip>
#include <sstream>
#include <string>

#include "arg_parse.hpp"

ArgParse::ArgParse(std::string desc)
{
    desc_ = desc;
}

std::string ArgParse::get_help()
{
    std::stringstream ss;

    ss << desc_ << " Help Page\n\n";
    ss << desc_ << '\n';
    for (auto& arg : args_)  // process required arguments
    {
        if (arg.second.required)
        {
            ss << "-" << arg.first << " ";
        }
    }
    bool first = true;
    for (auto& arg : args_)  // now process non-required arguments
    {
        if (!arg.second.required)
        {
            if (first)
            {
                ss << "[";
            }
            first = false;
            ss << "-" << arg.first << " ";
        }
    }
    if (!first)
    {
        ss << "]";
    }
    ss << '\n';
    ss << "Arguments:\n";
    for (auto& arg : args_)
    {
        ss << "\t" << arg.first << "/--" << arg.second.name << "\t" << arg.second.desc << '\n';
    }
    return(ss.str());
}

std::string ArgParse::get_param_str()
{
    const std::string DELIM = ".";
    std::stringstream ss;

    for (auto a = args_.begin(); a != args_.end(); a++)
    {
        std::cout << a->first << '\n';

        if (a->second.type != Type::STRING)
        {
            if (a != args_.begin())
            {
                ss << DELIM;
            }
            ss << a->first;

            switch (a->second.type)
            {
            case Type::FLAG:
                ss << get_flag(a->first);
                break;

            case Type::INT:
                ss << get_int(a->first);
                break;

            case Type::DOUBLE:
                ss << std::fixed << std::setprecision(2) << get_double(a->first);
                break;

            default:
                std::cerr << "Error: bad type\n";
                break;
            }
        }
    }
    return(ss.str());
}

void ArgParse::parse_args(int argc, char **argv)
{
    std::string opt_str;

    for (auto a = args_.begin(); a != args_.end(); a++)
    {
        opt_str.push_back(a->first);

        if (a->second.type != Type::FLAG)
        {
            opt_str.push_back(':');
        }
    }

    char o;

    while ((o = getopt(argc, argv, opt_str.c_str())) != -1)
    {
        if (args_.count(o) == 0)
        {
            std::cerr << "Error: unrecognized argument " << o << '\n';
            continue;
        }


        Arg&a = args_[o];
        a.set = true; // this flag is set in the command line
        switch (a.type)
        {
        case Type::FLAG:
            a.value = "1";
            break;

        case Type::INT:
            a.value = optarg;
            break;

        case Type::DOUBLE:
            a.value = optarg;
            break;

        case Type::STRING:
            a.value = optarg;
            break;

        default:
            std::cerr << "Error: bad type\n";
            break;
        }
    }

    // now check if all required arguments have been provided
    for (auto& arg : args_)
    {
        if (arg.second.required && !arg.second.set)
        {
            std::cerr << "======================================\n\n" << "Missing argument: " << arg.second.name << "\n\n" << "======================================\n\n";
            std::cerr << this->get_help() << '\n';
            exit(1);
        }
    }
}

bool ArgParse::add_flag(char c, std::string name, std::string desc = "", bool required = false)
{
    if (args_.count(c) > 0)
    {
        return(false);
    }

    Arg a = { Type::FLAG, std::move(name), std::move(desc), "0", required, false };
    args_[c] = a;

    return(true);
}

bool ArgParse::add_int(char c, std::string name,
                       int def, std::string desc = "", bool required = false)
{
    if (args_.count(c) > 0)
    {
        return(false);
    }

    Arg a = { Type::INT, std::move(name), std::move(desc), std::to_string(def), required, false };
    args_[c] = a;

    return(true);
}

bool ArgParse::add_double(char c, std::string name,
                          double def, std::string desc = "", bool required = false)
{
    if (args_.count(c) > 0)
    {
        return(false);
    }

    Arg a = { Type::DOUBLE, std::move(name), std::move(desc), std::to_string(def), required, false };
    args_[c] = a;

    return(true);
}

bool ArgParse::add_string(char c, std::string name,
                          std::string def, std::string desc = "", bool required = false)
{
    if (args_.count(c) > 0)
    {
        return(false);
    }

    Arg a = { Type::STRING, std::move(name), std::move(desc), def, required, false };
    args_[c] = a;

    return(true);
}

std::string ArgParse::get_name(char c)
{
    return(args_[c].name);
}

std::string ArgParse::get_desc(char c)
{
    return(args_[c].desc);
}

bool ArgParse::get_flag(char c)
{
    Arg&a = args_[c];

    return((a.value == "1") ? true : false);
}

int ArgParse::get_int(char c)
{
    Arg&a = args_[c];

    return(std::stoi(a.value));
}

double ArgParse::get_double(char c)
{
    Arg&a = args_[c];

    return(std::stod(a.value));
}

std::string ArgParse::get_string(char c)
{
    Arg&a = args_[c];

    return(a.value);
}

bool ArgParse::is_set(char c)
{
    Arg&a = args_[c];

    return(a.set);
}
