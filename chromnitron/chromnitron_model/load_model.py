import torch
import torch.nn as nn
import os

def init_model(config):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    from chromnitron_model.chromnitron_models import Chromnitron
    model = Chromnitron()
    model = model.to(device)
    return model

def load_model_weights(model, weights_path):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    state_dict_list = torch.load(weights_path, map_location = device)
    try:
        model.load_state_dict(state_dict_list)
    except:
        print(f'Ignore for LoRA weights: model weights do not match model architecture for {weights_path}')
        model.load_state_dict(state_dict_list, strict = False)
    return model

def load_chromnitron(config, cap):
    model_resource = config['model_resource']
    finetune_mode = config['inference_config']['inference']['use_finetune']
    base_model_weights_path = f'{model_resource["root"]}/{model_resource["base_weights"]}'
    per_cap_lora_weights_path = f'{model_resource["root"]}/{model_resource["per_cap_lora_weights"]}/{cap}.pt'
    weights_exist = os.path.exists(per_cap_lora_weights_path)

    model = init_model(config)
    model = load_model_weights(model, base_model_weights_path)
    assert finetune_mode in ['auto', 'enable', 'disable']
    if finetune_mode in ['auto', 'enable']:
        if weights_exist:
            print(f'Using finetuned LoRA model for {cap}')
            model = load_lora_pretrained(model, base_model_weights_path, lora_r = 4) # Load main weights
            model = load_model_weights(model, per_cap_lora_weights_path)
            print('LoRA weights loaded')
        else:
            if finetune_mode == 'enable':
                raise ValueError(f'Finetuned model for {cap} not found')
            else:
                print(f'Finetuned model for {cap} not found, using base model')
    model.eval()
    return model

def load_lora_pretrained(model, finetune_model_path, lora_r):
    import loralib as lora
    print(f'Loading model from {finetune_model_path} with LoRA')
    replace_layers_with_lora(model, lora_layer_factory, r = lora_r)
    lora.mark_only_lora_as_trainable(model)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    weights = torch.load(finetune_model_path, map_location = device)
    if 'model' in weights:
        weights = weights['model']
    load_state_dict_to_lora(model, weights)
    return model

def load_state_dict_to_lora(model, state_dict):
    model_state_dict = model.state_dict()
    model_key_conv_list = list(model_state_dict.keys())
    no_conv_to_conv_dict = {}
    for key in model_key_conv_list:
        if '.conv.' in key:
            no_conv_key = key.replace('.conv.', '.')
            no_conv_to_conv_dict[no_conv_key] = key
        else:
            no_conv_to_conv_dict[key] = key
    # Modify the state dict
    new_state_dict = {}
    for key in state_dict.keys():
        new_key = key
        new_state_dict[no_conv_to_conv_dict[new_key]] = state_dict[key]
    model.load_state_dict(new_state_dict, strict = False)

def replace_layers_with_lora(model, lora_layer_factory, r):
    for name, module in model.named_children():
        if len(list(module.children())) > 0:
            # Recursively apply to child modules
            replace_layers_with_lora(module, lora_layer_factory, r)
        else:
            # Replace the layer with a LoRA layer
            setattr(model, name, lora_layer_factory(module, r))
    return model

def lora_layer_factory(original_layer, r):
    import loralib as lora
    if isinstance(original_layer, nn.Linear):
        return lora.Linear(original_layer.in_features, original_layer.out_features, r = r)
    elif isinstance(original_layer, nn.Embedding):
        return lora.Embedding(original_layer.num_embeddings, original_layer.embedding_dim, r = r)
    elif isinstance(original_layer, nn.Conv1d):
        return lora.Conv1d(original_layer.in_channels, 
                           original_layer.out_channels, 
                           original_layer.kernel_size[0], 
                           r = r, 
                           stride=original_layer.stride[0], padding=original_layer.padding[0], dilation=original_layer.dilation[0])
    elif isinstance(original_layer, nn.ConvTranspose1d):
        return lora.ConvTranspose1d(original_layer.in_channels, 
                                    original_layer.out_channels, 
                                    original_layer.kernel_size[0], 
                                    r = r, 
                                    stride=original_layer.stride[0], padding=original_layer.padding[0], dilation=original_layer.dilation[0])
    else:
        return original_layer